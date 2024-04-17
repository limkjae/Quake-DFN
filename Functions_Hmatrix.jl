function GroupAndSort(Input_Segment_Part, HowManyDivision, 
                        ArrangePoint, InitialBlockIdx, BlockBegin, TotalHierarchyLevel, 
                        CurrentLevel, HierarchyLevelHistory)

    R = kmeans(Input_Segment_Part[:,1:3]', HowManyDivision; maxiter=200)#, display=:iter)
    a = assignments(R) # get the assignments of points to clusters
    c = counts(R) # get the cluster sizes
    M = R.centers # get the cluster centers

    GroupDistanceFromArrangePoint = norm.(eachcol(M .- ArrangePoint))
    GroupOrder = sortperm(GroupDistanceFromArrangePoint)

    ClusterCenter = M[:,GroupOrder]
    ClusterSize = c[GroupOrder]
    # GroupDistanceFromArrangePoint[GroupOrder]
    OriginalAssignment = collect(Int, 1:length(GroupDistanceFromArrangePoint))

    for OrderIndex = 1:length(GroupDistanceFromArrangePoint)
        OriginalAssignment[GroupOrder[OrderIndex]] = OrderIndex
    end
    
    Input_Segment_Part[:,20]= OriginalAssignment[a]  .+ (InitialBlockIdx - 1)
    Input_Segment_Part[:,21]= GroupDistanceFromArrangePoint[a]
    # Input_Segment=[Input_Segment OriginalAssignment[a] GroupDistanceFromArrangePoint[a]]

    Input_Segment_Part=sortslices(Input_Segment_Part, dims = 1, by = x -> x[20])

    AddedBlock_Range_Level = zeros(Int,HowManyDivision, 3 + TotalHierarchyLevel)
    AddedBlock_Ctr_Diam = zeros(HowManyDivision, 4)
    Diameter = zeros(1, HowManyDivision)
    for i = 1 : HowManyDivision
        if i==1
            Index_i =1
        else 
            Index_i = sum(ClusterSize[1:i-1]) + 1 
        end
        
        Index_f = sum(ClusterSize[1:i])
        Diameter[i] = sqrt((maximum(Input_Segment_Part[Index_i:Index_f,1]) - 
                        minimum(Input_Segment_Part[Index_i:Index_f,1]))^2 +        
                      (maximum(Input_Segment_Part[Index_i:Index_f,2]) - 
                        minimum(Input_Segment_Part[Index_i:Index_f,2]))^2 +            
                      (maximum(Input_Segment_Part[Index_i:Index_f,3]) - 
                        minimum(Input_Segment_Part[Index_i:Index_f,3]))^2)   

        AddedBlock_Range_Level[i,1] = Index_i + BlockBegin -1
        AddedBlock_Range_Level[i,2] = Index_f + BlockBegin -1
        AddedBlock_Range_Level[i,3] = CurrentLevel
        AddedBlock_Range_Level[i,4:end] = HierarchyLevelHistory
        AddedBlock_Range_Level[i,3 + CurrentLevel] = i
        
        AddedBlock_Ctr_Diam[i,1:3] = ClusterCenter[:,i]'
        AddedBlock_Ctr_Diam[i,4] = Diameter[i]
    end

    # println(AddedBlock_Range_Level)
    return Input_Segment_Part, ClusterCenter, GroupOrder, ClusterSize, Diameter, 
            AddedBlock_Range_Level, AddedBlock_Ctr_Diam
end


function GroupAndSort_AllLevel(HowManyDivisionEachLevel, TotalHierarchyLevel, MinimumElementsToCut,
                                ArrangePoint, FaultCount, Input_Segment)

Block_Range_Level = zeros(Int,1, 3 + TotalHierarchyLevel)
Block_Ctr_Diam = zeros(1, 4)
BlockCountPrevLevel = zeros(Int,TotalHierarchyLevel)
LevelBeign = 0
LevelEnd = 0
for CurrentLevel = 1:TotalHierarchyLevel
    print("Level ",CurrentLevel, "\n")
    if CurrentLevel == 1
        BlockBegin = 1
        BlockEnd =FaultCount
        InitialBlockIdx = Input_Segment[BlockBegin,20]
        Input_Segment[BlockEnd:end,20] = Input_Segment[BlockEnd:end,20] .+ (HowManyDivisionEachLevel -1 )

        HierarchyLevelHistory = Block_Range_Level[1,4:end]
        
        Input_Segment[BlockBegin:BlockEnd,:], ClusterCenter, GroupOrder, ClusterSize, Diameter, 
            AddedBlock_Range_Level, AddedBlock_Ctr_Diam = 
            GroupAndSort(Input_Segment[BlockBegin:BlockEnd,:], HowManyDivisionEachLevel, 
                            ArrangePoint, InitialBlockIdx, BlockBegin, TotalHierarchyLevel, 
                            CurrentLevel, HierarchyLevelHistory)

        Block_Range_Level = [Block_Range_Level; AddedBlock_Range_Level]
        Block_Ctr_Diam = [Block_Ctr_Diam; AddedBlock_Ctr_Diam]

        Block_Range_Level = Block_Range_Level[2:end,:]
        Block_Ctr_Diam = Block_Ctr_Diam[2:end,:]
        BlockCountPrevLevel[CurrentLevel] = length(Block_Ctr_Diam[:,1])
        LevelBeign =1 
        LevelEnd = length(Block_Ctr_Diam[:,1])
        
    else 
        for highLevelidx = LevelBeign : LevelEnd
        BlockBegin = Int(Block_Range_Level[highLevelidx,1])
        BlockEnd = Int(Block_Range_Level[highLevelidx,2])
        if BlockEnd - BlockBegin < MinimumElementsToCut
            DivisionThisTime = 1
        else
            DivisionThisTime = HowManyDivisionEachLevel
        end
            HierarchyLevelHistory = Block_Range_Level[highLevelidx,4:end]
            InitialBlockIdx = Input_Segment[BlockBegin,20]
            Input_Segment[BlockEnd:end,20] = Input_Segment[BlockEnd:end,20] .+ (DivisionThisTime -1 )

            Input_Segment[BlockBegin:BlockEnd,:], ClusterCenter, GroupOrder, ClusterSize, Diameter, 
                AddedBlock_Range_Level, AddedBlock_Ctr_Diam = 
                GroupAndSort(Input_Segment[BlockBegin:BlockEnd,:], DivisionThisTime, 
                                ArrangePoint, InitialBlockIdx, BlockBegin, TotalHierarchyLevel, CurrentLevel, HierarchyLevelHistory)

            Block_Range_Level = [Block_Range_Level; AddedBlock_Range_Level]
            Block_Ctr_Diam = [Block_Ctr_Diam; AddedBlock_Ctr_Diam]
        
        end
        BlockCountPrevLevel[CurrentLevel] = length(Block_Ctr_Diam[:,1]) - sum(BlockCountPrevLevel[1:CurrentLevel-1])
        LevelBeign =  sum(BlockCountPrevLevel[1:CurrentLevel-1]) +1 
        LevelEnd = length(Block_Ctr_Diam[:,1])
    end
end
  return Block_Ctr_Diam, Block_Range_Level, Input_Segment
end

############################## Build Hierarchy ####################################


function BuildHierarchy(Block_Range_Level, Block_Ctr_Diam, DistDiamRatioCrit, TotalHierarchyLevel) 

    Ranks =  zeros(Int,1)
    ElementRange_SR = zeros(Int, 1,4)
    for Level = 1:TotalHierarchyLevel
        Range_Level_Current = Block_Range_Level[findall(x-> x== Level, Block_Range_Level[:,3]),:]
        Ctr_Diam_Current = Block_Ctr_Diam[findall(x-> x== Level, Block_Range_Level[:,3]),:]
        BlockedElementCount = length(Range_Level_Current[:,1])
        BlockCount = BlockedElementCount^2
        DoneInHigherLevel = 0
    
        for i=1:BlockedElementCount
            for j=1:BlockedElementCount
                DoneInHigherLevel = 0
                for k in eachindex(ElementRange_SR[:,1])
                    if ElementRange_SR[k,1] <= Range_Level_Current[i,1] < ElementRange_SR[k,2] &&
                        ElementRange_SR[k,3] <= Range_Level_Current[j,1] < ElementRange_SR[k,4]
                        DoneInHigherLevel = 1
                        # println(ElementRange_SR)
                    end
                end
                if DoneInHigherLevel == 0
                    Distance = norm(Ctr_Diam_Current[i,1:3] - Ctr_Diam_Current[j,1:3])
                    Diameter = (Ctr_Diam_Current[i,4] + Ctr_Diam_Current[j,4])/2
                    if Distance/Diameter > DistDiamRatioCrit * 2
                        Ranks = [Ranks 10]
                        Added = [Range_Level_Current[i,1] Range_Level_Current[i,2] Range_Level_Current[j,1] Range_Level_Current[j,2]]
                        ElementRange_SR = [ElementRange_SR; Added] 
                    elseif Distance/Diameter > DistDiamRatioCrit * 1.5
                      Ranks = [Ranks 15]
                      Added = [Range_Level_Current[i,1] Range_Level_Current[i,2] Range_Level_Current[j,1] Range_Level_Current[j,2]]
                      ElementRange_SR = [ElementRange_SR; Added]                    

                    elseif Distance/Diameter > DistDiamRatioCrit 
                      Ranks = [Ranks 30]
                      Added = [Range_Level_Current[i,1] Range_Level_Current[i,2] Range_Level_Current[j,1] Range_Level_Current[j,2]]
                      ElementRange_SR = [ElementRange_SR; Added] 
                    
                    elseif Level == TotalHierarchyLevel
                        Ranks = [Ranks 0]
                        Added = [Range_Level_Current[i,1] Range_Level_Current[i,2] Range_Level_Current[j,1] Range_Level_Current[j,2]]
                        ElementRange_SR = [ElementRange_SR; Added] 
                    end
                end
            end
        end
    end
    Ranks = Ranks[2:end]
    ElementRange_SR = ElementRange_SR[2:end,:]
    return Ranks, ElementRange_SR

end

function HmatSolver(NetDisp, ShearStiffness_H, BlockCount, ElementRange_SR, FaultCount)
    Elastic_Load_Disp = zeros(FaultCount)
    # @sync begin
        # @distributed for Blockidx in 1:BlockCount  
        # @turbo for Blockidx in 1:BlockCount    
    
    for Blockidx in 1:BlockCount    
            Elastic_Load_Disp[ElementRange_SR[Blockidx,1]:ElementRange_SR[Blockidx,2]] = 
                Elastic_Load_Disp[ElementRange_SR[Blockidx,1]:ElementRange_SR[Blockidx,2]] + 
                ShearStiffness_H[Blockidx] * NetDisp[ElementRange_SR[Blockidx,3]:ElementRange_SR[Blockidx,4]]
        end
    
    return Elastic_Load_Disp
end

function HmatSolver_Pararllel(NetDisp, ShearStiffness_H, ElementRange_SR, FaultCount, Par_ElementDivision, ThreadCount)

    Elastic_Load_DispPart = zeros(FaultCount,ThreadCount)
    @sync begin

            # for i=1:ThreadCount; println(i);                
            # @time for Blockidx = Par_ElementDivision[i]+1 : Par_ElementDivision[i+1]
            
            Threads.@threads for i=1:ThreadCount 
             for Blockidx = Par_ElementDivision[i]+1 : Par_ElementDivision[i+1]
                Elastic_Load_DispPart[ElementRange_SR[Blockidx,1]:ElementRange_SR[Blockidx,2],i] = 
                Elastic_Load_DispPart[ElementRange_SR[Blockidx,1]:ElementRange_SR[Blockidx,2],i] + 
                    ShearStiffness_H[Blockidx] * NetDisp[ElementRange_SR[Blockidx,3]:ElementRange_SR[Blockidx,4]]
                    
            end
        end
    end
    # println("One Step Over")
   Elastic_Load_DispP = sum(Elastic_Load_DispPart, dims=2)
    return Elastic_Load_DispP
end


function ParallelOptimization(ShearStiffness_H, ElementRange_SR, FaultCount, Par_BlockCount, ThreadCount)
    println("Parallel Calculation Optimizing")
    RepeatCount=50
    NetDisp = ones(FaultCount)
    Elastic_Load_DispPart = zeros(FaultCount)

    # RecordedTime = ones(Par_BlockCount) * 1000
    RecordedTime = zeros(Par_BlockCount) 
    TimeEachRepeat = zeros(RepeatCount)
    # RecordedTime = zeros(Par_BlockCount) 
             for Blockidx = 1 : Par_BlockCount       
                for RepeatTime = 1:RepeatCount        
                    ElapsedTime= @elapsed begin
                    Elastic_Load_DispPart[ElementRange_SR[Blockidx,1]:ElementRange_SR[Blockidx,2]] = 
                    Elastic_Load_DispPart[ElementRange_SR[Blockidx,1]:ElementRange_SR[Blockidx,2]] + 
                        ShearStiffness_H[Blockidx] * NetDisp[ElementRange_SR[Blockidx,3]:ElementRange_SR[Blockidx,4]]
                    end
                    TimeEachRepeat[RepeatTime] = ElapsedTime
                    sort!(TimeEachRepeat)
                    # RecordedTime[Blockidx] = sum(TimeEachRepeat[round(Int,RepeatCount/10):round(Int,RepeatCount*2/10)])
                    RecordedTime[Blockidx] = sum(TimeEachRepeat[1:round(Int,RepeatCount*1/10)])
                    # if ElapsedTime < RecordedTime[Blockidx]
                    #     RecordedTime[Blockidx] = ElapsedTime
                    # end
                end
                    # RecordedTime[Blockidx] = RecordedTime[Blockidx] + ElapsedTime
            end 
                
        
    CumRecordTime = cumsum(RecordedTime)
    CumRecordTimeNorm = CumRecordTime / maximum(CumRecordTime)
    figure(3); plot(CumRecordTimeNorm)

    ElementGap = 1 / ThreadCount
    Par_ElementDivision = zeros(Int, ThreadCount +1)
    for i=1:ThreadCount
        if i==ThreadCount 
            Par_ElementDivision[i+1] = Par_BlockCount
        else 
            Par_ElementDivision[i+1] = findmin(abs.(CumRecordTimeNorm .- ElementGap * i))[2]
        end
    end
    println(Par_ElementDivision)
    # println(RecordedTime)
    # Elastic_Load_DispP = sum(Elastic_Load_DispPart[:,:], dims=2)
    return Par_ElementDivision
end

#= 
function GroupAndSort(Input_Segment_Part, HowManyDivision, 
    ArrangePoint, InitialBlockIdx, BlockBegin, TotalHierarchyLevel, CurrentLevel, HierarchyLevelHistory)



R = kmeans(Input_Segment_Part[:,1:3]', HowManyDivision; maxiter=200, display=:iter)

a = assignments(R) # get the assignments of points to clusters
c = counts(R) # get the cluster sizes
M = R.centers # get the cluster centers


GroupDistanceFromArrangePoint = norm.(eachcol(M .- ArrangePoint))
GroupOrder = sortperm(GroupDistanceFromArrangePoint)

ClusterCenter = M[:,GroupOrder]
ClusterSize = c[GroupOrder]
# GroupDistanceFromArrangePoint[GroupOrder]
OriginalAssignment = collect(Int, 1:length(GroupDistanceFromArrangePoint))

for OrderIndex = 1:length(GroupDistanceFromArrangePoint)
OriginalAssignment[GroupOrder[OrderIndex]] = OrderIndex
end

Input_Segment_Part[:,20]= OriginalAssignment[a]  .+ (InitialBlockIdx - 1)
Input_Segment_Part[:,21]= GroupDistanceFromArrangePoint[a]
# Input_Segment=[Input_Segment OriginalAssignment[a] GroupDistanceFromArrangePoint[a]]

Input_Segment_Part=sortslices(Input_Segment_Part, dims = 1, by = x -> x[20])

AddedBlockInfo = zeros(HowManyDivision,6 + TotalHierarchyLevel)
AddedBlock_Range_Level = zeros(Int,HowManyDivision, 2 + TotalHierarchyLevel)
AddedBlock_Ctr_Diam = zeros(HowManyDivision, 4)
Diameter = zeros(1, 4)
for i = 1:HowManyDivision
if i==1
Index_i =1
else 
Index_i = sum(ClusterSize[1:i-1]) + 1 
end

Index_f = sum(ClusterSize[1:i])
Diameter[i] = sqrt((maximum(Input_Segment_Part[Index_i:Index_f,1]) - 
    minimum(Input_Segment_Part[Index_i:Index_f,1]))^2 +        
  (maximum(Input_Segment_Part[Index_i:Index_f,2]) - 
    minimum(Input_Segment_Part[Index_i:Index_f,2]))^2 +            
  (maximum(Input_Segment_Part[Index_i:Index_f,3]) - 
    minimum(Input_Segment_Part[Index_i:Index_f,3]))^2)   

AddedBlock_Range_Level[i,1] = Index_i + BlockBegin -1
AddedBlock_Range_Level[i,2] = Index_f + BlockBegin -1
AddedBlock_Range_Level[i,3:end] = HierarchyLevelHistory
AddedBlock_Range_Level[i,2 + CurrentLevel] = i

AddedBlock_Ctr_Diam[i,1:3] = ClusterCenter[:,i]'
AddedBlock_Ctr_Diam[i,4] = Diameter[i]
end


return Input_Segment_Part, ClusterCenter, GroupOrder, ClusterSize, Diameter, 
AddedBlock_Range_Level, AddedBlock_Ctr_Diam
end
=#