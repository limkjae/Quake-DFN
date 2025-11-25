function GroupAndSort(Input_Segment_Part,  HowManyDivision, 
                        ArrangePoint, InitialBlockIdx, BlockBegin, TotalHierarchyLevel, 
                        CurrentLevel, HierarchyLevelHistory, RorT)
    if RorT == "R"
        FaultCenter = Input_Segment_Part[:,1:3]'
        LevelIdx = 20
    else
        FaultCenterX = (Input_Segment_Part[:,1]+Input_Segment_Part[:,4]+Input_Segment_Part[:,7]) ./ 3
        FaultCenterY = (Input_Segment_Part[:,2]+Input_Segment_Part[:,5]+Input_Segment_Part[:,8]) ./ 3
        FaultCenterZ = (Input_Segment_Part[:,3]+Input_Segment_Part[:,6]+Input_Segment_Part[:,9]) ./ 3
        FaultCenter = [FaultCenterX FaultCenterY FaultCenterZ]'
        LevelIdx = 21
    end

    R = kmeans(FaultCenter, HowManyDivision; maxiter=200)#, display=:iter)
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
    
    Input_Segment_Part[:,LevelIdx]= OriginalAssignment[a]  .+ (InitialBlockIdx - 1)
    Input_Segment_Part[:,LevelIdx+1]= GroupDistanceFromArrangePoint[a]
    # Input_Segment=[Input_Segment OriginalAssignment[a] GroupDistanceFromArrangePoint[a]]

    Input_Segment_Part=sortslices(Input_Segment_Part, dims = 1, by = x -> x[LevelIdx])

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


function GroupAndSort_AllLevel(TotalHierarchyLevel, MinimumElementsToCut,
                                ArrangePoint, Input_Segment,RorT)
if RorT == "R"
    VconstIdx = 16; LevelIdx = 20
else
    VconstIdx = 19; LevelIdx = 21
end

Block_Range_Level = zeros(Int,1, 3 + TotalHierarchyLevel)
Block_Ctr_Diam = zeros(1, 4)
BlockCountPrevLevel = zeros(Int,TotalHierarchyLevel)
LevelBeign = 0
LevelEnd = 0
LoadingFaultExist = 0
NonLoadingFaultRange = [minimum(findall(x->x==0, Input_Segment[:,VconstIdx])), maximum(findall(x->x==0, Input_Segment[:,VconstIdx]))]
if maximum(Input_Segment[:,1VconstIdx])>0
    LoadingFaultExist= 1
    LoadingFaultRange = [minimum(findall(x->x>0, Input_Segment[:,VconstIdx])), maximum(findall(x->x>0, Input_Segment[:,VconstIdx]))]
    Block_Range_Level[1:2] = LoadingFaultRange # Loading fault range level
    Block_Range_Level[3] = 0
end

for CurrentLevel = 1:TotalHierarchyLevel
    print("Level ",CurrentLevel, "\n")
    if CurrentLevel == 1
        BlockBegin = NonLoadingFaultRange[1]
        BlockEnd = NonLoadingFaultRange[2]
        InitialBlockIdx = Input_Segment[BlockBegin,LevelIdx]
        Input_Segment[BlockEnd:end,LevelIdx] = Input_Segment[BlockEnd:end,LevelIdx] .+ (2 -1 )

        HierarchyLevelHistory = Block_Range_Level[1,4:end]
        
        Input_Segment[BlockBegin:BlockEnd,:], ClusterCenter, GroupOrder, ClusterSize, Diameter, 
            AddedBlock_Range_Level, AddedBlock_Ctr_Diam = 
            GroupAndSort(Input_Segment[BlockBegin:BlockEnd,:], 2,
                            ArrangePoint, InitialBlockIdx, BlockBegin, TotalHierarchyLevel, 
                            CurrentLevel, HierarchyLevelHistory, RorT)

        Block_Range_Level = [Block_Range_Level; AddedBlock_Range_Level]
        Block_Ctr_Diam = [Block_Ctr_Diam; AddedBlock_Ctr_Diam]
        BlockCountPrevLevel[CurrentLevel] = length(Block_Ctr_Diam[:,1])

        # if LoadingFaultExist == 0 
        #     Block_Range_Level = Block_Range_Level[2:end,:]
        #     Block_Ctr_Diam = Block_Ctr_Diam[2:end,:]
        #     LevelBeign =1 
        #     LevelEnd = length(Block_Ctr_Diam[:,1])
        # else 
            LevelBeign =2
            LevelEnd = length(Block_Ctr_Diam[:,1])
        # end
        
    else 
        for highLevelidx = LevelBeign : LevelEnd
            BlockBegin = Int(Block_Range_Level[highLevelidx,1])
            BlockEnd = Int(Block_Range_Level[highLevelidx,2])
            if BlockEnd - BlockBegin < MinimumElementsToCut
                DivisionThisTime = 1
            else
                DivisionThisTime = 2
            end
                HierarchyLevelHistory = Block_Range_Level[highLevelidx,4:end]
                InitialBlockIdx = Input_Segment[BlockBegin,LevelIdx]
                Input_Segment[BlockEnd:end,LevelIdx] = Input_Segment[BlockEnd:end,LevelIdx] .+ (DivisionThisTime -1 )

                Input_Segment[BlockBegin:BlockEnd,:], ClusterCenter, GroupOrder, ClusterSize, Diameter, 
                    AddedBlock_Range_Level, AddedBlock_Ctr_Diam = 
                    GroupAndSort(Input_Segment[BlockBegin:BlockEnd,:], DivisionThisTime, 
                                    ArrangePoint, InitialBlockIdx, BlockBegin, TotalHierarchyLevel, CurrentLevel, HierarchyLevelHistory, RorT) 

                Block_Range_Level = [Block_Range_Level; AddedBlock_Range_Level]
                Block_Ctr_Diam = [Block_Ctr_Diam; AddedBlock_Ctr_Diam]
            
        end
        BlockCountPrevLevel[CurrentLevel] = length(Block_Ctr_Diam[:,1]) - sum(BlockCountPrevLevel[1:CurrentLevel-1])
        # if LoadingFaultExist == 0 
        #     LevelBeign =  sum(BlockCountPrevLevel[1:CurrentLevel-1]) 
        # else
        #     LevelBeign =  sum(BlockCountPrevLevel[1:CurrentLevel-1]) +1 
        # end
        LevelBeign =  sum(BlockCountPrevLevel[1:CurrentLevel-1]) +1 
        LevelEnd = length(Block_Ctr_Diam[:,1])
    end
end
  return Block_Ctr_Diam, Block_Range_Level, Input_Segment, LoadingFaultExist
end

############################## Build Hierarchy ####################################


function BuildHierarchy(Block_Range_Level, Block_Ctr_Diam, DistDiamRatioCrit, TotalHierarchyLevel, LoadingFaultExist) 

    Admissible =  zeros(Int,1)
    ElementRange_SR = zeros(Int, 1,4)

    if LoadingFaultExist ==1 
        Admissible = [Admissible 0]
        Added = [Block_Range_Level[1,1]     Block_Range_Level[1,2]     1     Block_Range_Level[1,1]-1 ]
        ElementRange_SR = [ElementRange_SR; Added] 
    end

    for Level = 1:TotalHierarchyLevel
        Range_Level_Current = Block_Range_Level[findall(x-> x== Level, Block_Range_Level[:,3]),:]
        Ctr_Diam_Current = Block_Ctr_Diam[findall(x-> x== Level, Block_Range_Level[:,3]),:]
        BlockedElementCount = length(Range_Level_Current[:,1])
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
                    Diameter = maximum([Ctr_Diam_Current[i,4], Ctr_Diam_Current[j,4]])
                    # Diameter = (Ctr_Diam_Current[i,4] + Ctr_Diam_Current[j,4])/2
                    if Distance/Diameter > DistDiamRatioCrit 
                        Admissible = [Admissible 1]
                        Added = [Range_Level_Current[i,1] Range_Level_Current[i,2] Range_Level_Current[j,1] Range_Level_Current[j,2]]
                        ElementRange_SR = [ElementRange_SR; Added] 
                   
                    elseif Level == TotalHierarchyLevel
                        Admissible = [Admissible 0]
                        Added = [Range_Level_Current[i,1] Range_Level_Current[i,2] Range_Level_Current[j,1] Range_Level_Current[j,2]]
                        ElementRange_SR = [ElementRange_SR; Added] 
                    end
                end
            end
        end
    end
    Admissible = Admissible[2:end]
    ElementRange_SR = ElementRange_SR[2:end,:]
    return Admissible, ElementRange_SR

end

function HmatSolver(NetDisp, ShearStiffness_H, BlockCount, ElementRange_SR, FaultCount)
    Elastic_Load_Disp = zeros(FaultCount)
  
        for Blockidx in 1:BlockCount    
            Elastic_Load_Disp[ElementRange_SR[Blockidx,1]:ElementRange_SR[Blockidx,2]] = 
                Elastic_Load_Disp[ElementRange_SR[Blockidx,1]:ElementRange_SR[Blockidx,2]] + 
                ShearStiffness_H[Blockidx] * NetDisp[ElementRange_SR[Blockidx,3]:ElementRange_SR[Blockidx,4]]
        end
    return Elastic_Load_Disp
end






function HmatSolver_Pararllel_Old(NetDisp, ShearStiffness_H, ElementRange_SR, FaultCount, Par_ElementDivision, ThreadCount)

    Elastic_Load_DispPart = zeros(FaultCount,ThreadCount)    
    
    @sync begin
        Threads.@threads for i=1:ThreadCount 
            Elastic_Load_DispPart[:,i] = SolveEachThread(Par_ElementDivision[i]+1, Par_ElementDivision[i+1], ElementRange_SR, ShearStiffness_H, NetDisp, FaultCount)
        end
    end
    Elastic_Load_DispP = sum(Elastic_Load_DispPart, dims=2)
  
    return Elastic_Load_DispP
end




function HmatSolver_ThreadTime(NetDisp, ShearStiffness_H, ElementRange_SR, FaultCount, Par_ElementDivision, ThreadCount)

    Elastic_Load_DispPart = zeros(FaultCount,ThreadCount)
    
        for i=1:ThreadCount; #println(i);              
            @time Elastic_Load_DispPart[:,i] = SolveEachThread(Par_ElementDivision[i]+1, Par_ElementDivision[i+1], ElementRange_SR, ShearStiffness_H, NetDisp, FaultCount)
        end  

    # @sync begin
    #     Threads.@threads for i=1:ThreadCount 
    #          Elastic_Load_DispPart[:,i] = SolveEachThread(Par_ElementDivision[i]+1, Par_ElementDivision[i+1], ElementRange_SR, ShearStiffness_H, NetDisp, FaultCount)
    #     end
    # end


    Elastic_Load_DispP = sum(Elastic_Load_DispPart, dims=2)
    return Elastic_Load_DispP
end





function HmatSolver_SpeedTest(NetDisp, ShearStiffness_H, ElementRange_SR, FaultCount, Par_ElementDivision, ThreadCount)
    TestRep = 5
    Elastic_Load_DispPart = zeros(FaultCount,ThreadCount)
    ElapsedTime= zeros(ThreadCount,TestRep)
    ElapsedTMin= zeros(ThreadCount)
    for TestIndex = 1:TestRep
        for i=1:ThreadCount; #println(i);              
            ElapsedTime[i,TestIndex] = @elapsed begin 
                Elastic_Load_DispPart[:,i] = SolveEachThread(Par_ElementDivision[i]+1, Par_ElementDivision[i+1], 
                                                ElementRange_SR, 
                                                ShearStiffness_H,
                                                NetDisp, FaultCount)
            end
        end  
    end
    ElapsedTMin = vec(median(ElapsedTime, dims=2))

    return ElapsedTMin
end



function SolveEachThread(BlockI, BlockF, ElementRange_SR, ShearStiffness_H, NetDisp, FaultCount)
    Elastic_Load_D = zeros(FaultCount)
    for Blockidx = BlockI:BlockF
        Elastic_Load_D[ElementRange_SR[Blockidx,3]:ElementRange_SR[Blockidx,4]] +=
         ShearStiffness_H[Blockidx] * (NetDisp[ElementRange_SR[Blockidx,1]:ElementRange_SR[Blockidx,2]])
    end
    return Elastic_Load_D
end


function  ParallelOptimization(ShearStiffness_H, ElementRange_SR, 
    FaultCount, BlockCount, ThreadCount, MaxRatioAllowed, MaxIteration)
    println("Parallel Calculation Optimizing")
    ############## calculating initial guess for division ##############
    RepeatCount=10
    NetDisp = ones(FaultCount)
    Elastic_Load_DispPart = zeros(FaultCount)

    RecordedTime = zeros(BlockCount) 
    TimeEachRepeat = zeros(RepeatCount)
        for Blockidx = 1 : BlockCount       
            for RepeatTime = 1:RepeatCount        
                ElapsedTime= @elapsed begin
                    Elastic_Load_DispPart[ElementRange_SR[Blockidx,3]:ElementRange_SR[Blockidx,4]] += 
                    ShearStiffness_H[Blockidx] * NetDisp[ElementRange_SR[Blockidx,1]:ElementRange_SR[Blockidx,2]]
                end
                TimeEachRepeat[RepeatTime] = ElapsedTime
            end
            sort!(TimeEachRepeat)
            RecordedTime[Blockidx] = sum(TimeEachRepeat[1:round(Int,RepeatCount*3/10)])
        end 
                
        
    CumRecordTime = cumsum(RecordedTime)
    CumRecordTimeNorm = CumRecordTime / maximum(CumRecordTime)


    Par_ElementDivision = zeros(Int, ThreadCount +1)
    for i=1:ThreadCount
        if i==ThreadCount 
            Par_ElementDivision[i+1] = BlockCount
        else 
            Par_ElementDivision[i+1] = findmin(abs.(CumRecordTimeNorm .- 1 / ThreadCount * i))[2]
        end
    end
    println("First Guess: \n", Par_ElementDivision)
  

    ################## Find the best division by actual computations and adjustment ################

    ElementCountPerDivision = zeros(Int,ThreadCount)
    for i=1:ThreadCount
    ElementCountPerDivision[i] = Par_ElementDivision[i+1] - Par_ElementDivision[i] 
    end
    
    
    ElapsedTime = HmatSolver_SpeedTest(NetDisp, ShearStiffness_H, ElementRange_SR, FaultCount, Par_ElementDivision, ThreadCount)
    AverageTime = mean(ElapsedTime)
    
    Minimum_Par_ElementDivision = zeros(length(Par_ElementDivision))
    Minimum_ElapsedTimeRatio = 1e5
    iteration=0
    Termination = 1

    while Termination ==1
        iteration=iteration+1
        for i=1:ThreadCount
            if ElapsedTime[i] > AverageTime
                ElapsedOverAverage = ElapsedTime[i] / AverageTime
                ElementCountPerDivision[i] = round(Int,ElementCountPerDivision[i] * (1 - (1- 1/ElapsedOverAverage)*0.5))
            else 
                ElapsedOverAverage = ElapsedTime[i] / AverageTime
                if ElementCountPerDivision[i] ==1 
                    ElementCountPerDivision[i] = round(Int,ElementCountPerDivision[i] * (1 - (1- 1/ElapsedOverAverage)))
                else
                ElementCountPerDivision[i] = round(Int,ElementCountPerDivision[i] * (1 - (1- 1/ElapsedOverAverage)*0.03))
                end
            end        
            if ElementCountPerDivision[i] < 1
                ElementCountPerDivision[i] = 1
            end
        end
    
        ElementCountPerDivision = round.(Int, ElementCountPerDivision * BlockCount / sum(ElementCountPerDivision))
    
        for i=1:ThreadCount -1 
            Par_ElementDivision[i+1] = ElementCountPerDivision[i] + Par_ElementDivision[i]
        end
        Par_ElementDivision[end]= BlockCount

        for i=1:ThreadCount -1 
            if Par_ElementDivision[end-i] > Par_ElementDivision[end-i+1] 
                Par_ElementDivision[end-i] = Par_ElementDivision[end-i+1] 
            end
        end


        for i=1:ThreadCount
            ElementCountPerDivision[i] = Par_ElementDivision[i+1] - Par_ElementDivision[i] 
        end
    
    
        ElapsedTime = HmatSolver_SpeedTest(NetDisp, ShearStiffness_H, ElementRange_SR, FaultCount, Par_ElementDivision, ThreadCount)
    
        Maxindex = argmax(ElapsedTime, dims=1)[1]
        AverageTime = mean(ElapsedTime)
        maxAverageTimeRatio = maximum(ElapsedTime/AverageTime)
        println(iteration, ", Maxtime / Average: ", maxAverageTimeRatio)
    
        if maxAverageTimeRatio < Minimum_ElapsedTimeRatio
            Minimum_ElapsedTimeRatio = copy(maxAverageTimeRatio)
            Minimum_Par_ElementDivision = copy(Par_ElementDivision)
        end
    
        if maxAverageTimeRatio < MaxRatioAllowed
            Termination = 0
            println(Par_ElementDivision)
        elseif iteration == MaxIteration
            Termination = 0
            Par_ElementDivision = copy(Minimum_Par_ElementDivision)
            println(Par_ElementDivision)
        end
    
        # figure(11)
        # plot(ElapsedTime)
        # plot([0,ThreadCount],[AverageTime,AverageTime])
        # figure(12)
        # plot(log10.(ElementCountPerDivision))
            
    end
    return Par_ElementDivision
end


function  ParallelOptimization2(Stiffness_H, BlockCount, ThreadCount, Ranks, Admissible)
    println("Parallel Calculation Optimizing")
    
    Complexity = zeros(BlockCount)
    for BlockIdx=1:BlockCount
        if Admissible[BlockIdx] > 0
            if Ranks[BlockIdx] == 0; Ranks[BlockIdx] = 1; end
            Complexity[BlockIdx] = (size(Stiffness_H[BlockIdx])[1] + size(Stiffness_H[BlockIdx])[2]) * Ranks[BlockIdx]
        else
            Complexity[BlockIdx] = (size(Stiffness_H[BlockIdx])[1] * size(Stiffness_H[BlockIdx])[2]) 
        end        
    end

    CumSum = cumsum(Complexity) .- 1e5
    Division = CumSum[end] / ThreadCount
    Par_ElementDivision = [0]
    CurrentDiv = 1
    for BlockIdx=1:BlockCount

        if CumSum[BlockIdx] > Division * CurrentDiv
            Par_ElementDivision = [Par_ElementDivision; BlockIdx]
            CurrentDiv += 1
        end 

    end

    Par_ElementDivision = [Par_ElementDivision; BlockCount]


    return Par_ElementDivision

end

function StiffnessTransitionToLoading(ShearStiffness_H, ElementRange_SR, FaultCount)
    TotalBlock = size(ElementRange_SR, 1)
    K_Self = zeros(FaultCount)
    LoadingSiffnessH = Array{Any}(undef, TotalBlock)
    for i=1:TotalBlock 
        LoadingSiffnessH[i] = copy(ShearStiffness_H[i])
        if ElementRange_SR[i,1] == ElementRange_SR[i,3] && 
            ElementRange_SR[i,2] == ElementRange_SR[i,4]
            ElementStartsAt = ElementRange_SR[i,1]
            for j = 1:size(ShearStiffness_H[i],1)
                K_Self[ElementStartsAt+j-1] = -ShearStiffness_H[i][j,j] 
                
                LoadingSiffnessH[i][j,j] = 0      
                
            end
        end
    end
    return LoadingSiffnessH, K_Self
end


function GetKself(ShearStiffness_H, ElementRange_SR, FaultCount)
    TotalBlock = size(ElementRange_SR, 1)
    K_Self = zeros(FaultCount)
    for i=1:TotalBlock 
        if ElementRange_SR[i,1] == ElementRange_SR[i,3] && 
            ElementRange_SR[i,2] == ElementRange_SR[i,4]
            ElementStartsAt = ElementRange_SR[i,1]
            for j = 1:size(ShearStiffness_H[i],1)
                K_Self[ElementStartsAt+j-1] = -ShearStiffness_H[i][j,j] 
            end
        end
    end
    return K_Self
end


function HMatBlockPlot(FaultCenter,FaultLengthStrike, FaultLengthDip, FaultStrikeAngle,
    FaultDipAngle, FaultLLRR, InputProperty, PlotRotation, MinMax_Axis, ColorMinMax, Transparent, Edge, LoadingFaultCount)

    # FaultPlot_3D_Color_General(FaultCenter[1:FaultCount-LoadingFaultCount,:],FaultLengthStrike[1:FaultCount-LoadingFaultCount], FaultLengthDip[1:FaultCount-LoadingFaultCount],
    # FaultStrikeAngle[1:FaultCount-LoadingFaultCount], FaultDipAngle[1:FaultCount-LoadingFaultCount], FaultLLRR[1:FaultCount-LoadingFaultCount],ReMidPressure[:,1:FaultCount-LoadingFaultCount], 
    # PlotRotation, MinMax_Axis, Transparent)

    cm = get_cmap(:tab10)
    
    buffer=300
    xmin=minimum(FaultCenter[1:end-LoadingFaultCount,1])-buffer
    xmax=maximum(FaultCenter[1:end-LoadingFaultCount,1])+buffer
    ymin=minimum(FaultCenter[1:end-LoadingFaultCount,2])-buffer
    ymax=maximum(FaultCenter[1:end-LoadingFaultCount,2])+buffer
    zmax=0
    zmin=-maximum(FaultCenter[1:end-LoadingFaultCount,3])-buffer

    xCenter=(xmin+xmax)/2
    yCenter=(ymin+ymax)/2
    zCenter=(zmin+zmax)/2
    
    maxlength=maximum([xmax-xmin,ymax-ymin, -zmin+zmax])

    if MinMax_Axis == 0
        xmin=xCenter-maxlength/2
        xmax=xCenter+maxlength/2
        ymin=yCenter-maxlength/2
        ymax=yCenter+maxlength/2
        zmin=zCenter-maxlength/2

    else
        xmin=MinMax_Axis[1,1]
        xmax=MinMax_Axis[1,2]
        ymin=MinMax_Axis[2,1]
        ymax=MinMax_Axis[2,2]
        zmin=MinMax_Axis[3,1]
        zmax=MinMax_Axis[3,2]
    end
    if ColorMinMax==0
    MaxValue=maximum(InputProperty)
    MinValue=minimum(InputProperty)
        if MaxValue == MinValue
            ColorRange = MaxValue*0.1
            MaxValue = MaxValue + ColorRange
            MinValue = MinValue - ColorRange
        end
    else
    MaxValue=ColorMinMax[2]
    MinValue=ColorMinMax[1]
    end
    #MaxDisp=maximum(ReMidDisp)
    #ReMidDisp_Plot=copy(ReMidDisp)
    # for FaultIdx in eachindex(FaultLengthStrike)  
    for FaultIdx = 1: length(FaultLengthStrike)  - LoadingFaultCount
    # for FaultIdx = 1125
        #print(FaultIdx)
        RotMatStrike=[cosd(FaultStrikeAngle[FaultIdx]) -sind(FaultStrikeAngle[FaultIdx]) 0
            sind(FaultStrikeAngle[FaultIdx]) cosd(FaultStrikeAngle[FaultIdx]) 0
            0  0 1]
        RotMatDip=[1 0  0
        0 cosd(FaultDipAngle[FaultIdx]) -sind(FaultDipAngle[FaultIdx])
        0 sind(FaultDipAngle[FaultIdx]) cosd(FaultDipAngle[FaultIdx])]
                
        p1=RotMatStrike*RotMatDip*[FaultLengthStrike[FaultIdx]/2;-FaultLengthDip[FaultIdx]/2;0] + [FaultCenter[FaultIdx,1]; FaultCenter[FaultIdx,2]; -FaultCenter[FaultIdx,3]];
        p2=RotMatStrike*RotMatDip*[-FaultLengthStrike[FaultIdx]/2;-FaultLengthDip[FaultIdx]/2;0] + [FaultCenter[FaultIdx,1]; FaultCenter[FaultIdx,2]; -FaultCenter[FaultIdx,3]];
        p3=RotMatStrike*RotMatDip*[-FaultLengthStrike[FaultIdx]/2;+FaultLengthDip[FaultIdx]/2;0]+ [FaultCenter[FaultIdx,1]; FaultCenter[FaultIdx,2]; -FaultCenter[FaultIdx,3]];
        p4=RotMatStrike*RotMatDip*[FaultLengthStrike[FaultIdx]/2;+FaultLengthDip[FaultIdx]/2;0]+ [FaultCenter[FaultIdx,1]; FaultCenter[FaultIdx,2]; -FaultCenter[FaultIdx,3]];
        
        fig = figure(1)
        # PlotValue=(InputProperty[FaultIdx]-MinValue)/(MaxValue-MinValue)
        PlotValue=(InputProperty[FaultIdx] % 10)/10
        #if FaultIdx==length(FaultLengthStrike)
        #    ion()
        #end
        #ion() # ioff() for turn off drawnow
        art3d = PyObject(PyPlot.art3D)

        verts2 = ([tuple(p1...); tuple(p2...); tuple(p3...); tuple(p4...)],)
        #p3c = PyObject(art3d.Poly3DCollection(verts2, linewidths=1, alpha=ReMidVel[PlotStep,FaultIdx]/MaxVel))
        #p3c = PyObject(art3d.Poly3DCollection(verts2, linewidths=1, alpha=0.05+ReMidVel[PlotStep,FaultIdx]/MaxVel*0.95)) 
        #p3c = PyObject(art3d.Poly3DCollection(verts2, linewidths=1, alpha=0.05+ReMidVelLog[PlotStep,FaultIdx]/MaxVelLog*0.95))
        p3c = PyObject(art3d.Poly3DCollection(verts2, linewidths=1))
        # ReMidVel[PlotStep,FaultIdx]/MaxVel))

        ax = gca(projection="3d")
        pycall(ax.add_collection3d, PyAny, p3c)
        xlim(xmin,xmax )
        ylim(ymin,ymax )
        zlim(zmin,zmax )

        if Transparent == 0
            face_color = [cm(PlotValue)[1], cm(PlotValue)[2],cm(PlotValue)[3],1]#PlotValue]#ReMidVel[PlotStep,FaultIdx]]
            edge_color = [0.5, 0.5, 0.5, 1.0]#ReMidVel[PlotStep,FaultIdx]]
        else 
            face_color = [cm(PlotValue)[1], cm(PlotValue)[2],cm(PlotValue)[3],0.3]#PlotValue]#ReMidVel[PlotStep,FaultIdx]]
            edge_color = [0.5, 0.5, 0.5, 1.0]#ReMidVel[PlotStep,FaultIdx]]

        end
        if Edge == 0
        edge_color = [0 0 0 0]
        end
        
        pycall(p3c.set_facecolor, PyAny, face_color)
        pycall(p3c.set_edgecolor, PyAny, edge_color)
        ax.view_init(PlotRotation[1],PlotRotation[2])

    end
    return MaxValue, MinValue


end

function SolveAx_b(LoadingStiffnessH, K_Self, InitialShearStress, ElementRange_SR, FaultCount, 
                Par_ElementDivision, ThreadCount, Epsilon_MaxDiffRatio)
    ########### Find Initial Displacement by Gauss Seidel ############
    println("Iteratively Calculating Initial Displacement")
    DispI_k = ones(FaultCount)
    MaxDiff = 1
    iteration=0
    while MaxDiff > Epsilon_MaxDiffRatio
        iteration = iteration+1
        if iteration % 20 ==0
        println("Iteration until 1 > ",MaxDiff / Epsilon_MaxDiffRatio)
        end
    EFTerm = HmatSolver_Pararllel(DispI_k, LoadingStiffnessH, ElementRange_SR, Par_ElementDivision, ThreadCount, zeros(FaultCount, ThreadCount))
    
    DispI_kp1 = (EFTerm + InitialShearStress) ./ K_Self
    MaxDiff = maximum(abs.((DispI_k - DispI_kp1) ./ DispI_k))
    DispI_k=copy(DispI_kp1)
    end
    return DispI_k
end



function HmatSolver_Pararllel(NetDisp, Stiffness_H, ElementRange_SR, 
    Par_ElementDivision, ThreadCount, Elastic_Load_EachThread)
            
        @inbounds Threads.@threads for ThreadIdx in 1:ThreadCount 

            for Blockidx = Par_ElementDivision[ThreadIdx]+1:Par_ElementDivision[ThreadIdx+1]
                @inbounds @views Elastic_Load_EachThread[ElementRange_SR[Blockidx,3]:ElementRange_SR[Blockidx,4], ThreadIdx] +=
                    Stiffness_H[Blockidx] * NetDisp[ElementRange_SR[Blockidx,1]:ElementRange_SR[Blockidx,2]]
            end   
      
                
        end
        Elastic_Load_DispP = sum(Elastic_Load_EachThread, dims=2)

    return Elastic_Load_DispP
end

function SolveAx_b_Jacobi(LoadingStiffnessH, K_Self, InitialShearStress, ElementRange_SR, FaultCount, 
                            Par_ElementDivision, ThreadCount, Epsilon_MaxDiffRatio)
    ########### Find Initial Displacement by Gauss Seidel ############
    println("Iteratively Calculating Initial Displacement")
    DispI_k = ones(FaultCount)
    MaxDiff = 1
    iteration=0
    println("Jabobi")
    while MaxDiff > Epsilon_MaxDiffRatio
    iteration = iteration+1
    if iteration % 20 ==0
    println("Iteration until 1 > ",MaxDiff / Epsilon_MaxDiffRatio)
    end
    
    EFTerm = HmatSolver_Pararllel(DispI_k, LoadingStiffnessH, ElementRange_SR, Par_ElementDivision, ThreadCount, zeros(FaultCount, ThreadCount))

    DispI_kp1 = (EFTerm + InitialShearStress) ./ K_Self
    MaxDiff = maximum(abs.((DispI_k - DispI_kp1)))
    DispI_k=copy(DispI_kp1)
    end
    return DispI_k
end



function SolveAx_b_GaussSeidel(LoadingStiffnessH, K_Self, InitialShearStress, ElementRange_SR, 
    FaultCount, Par_ElementDivision_Shear, ThreadCount, DispCritEps, Ranks_Shear, w_factor)


    BlockCount = length(Ranks_Shear)
    UpperBlockCount = 0
    LowerBlockCount = 0
    FullBlockCount = 0
    Stiffness_L = Any[]
    Stiffness_U =Any[]
    ElementRange_L = zeros(Int,1,4)
    ElementRange_U = zeros(Int,1,4)
    Rank_L = zeros(Int,1)
    Rank_U = zeros(Int,1)

    BlockSolution = Any[]

    for BlockIndex = 1:BlockCount
        BlockElementCount = ElementRange_SR[BlockIndex,4] -ElementRange_SR[BlockIndex,3] +1
        if ElementRange_SR[BlockIndex,1] == ElementRange_SR[BlockIndex,3] && ElementRange_SR[BlockIndex,2] == ElementRange_SR[BlockIndex,4]

            UpperBlockCount += 1
            UpperMatrix = copy(LoadingStiffnessH[BlockIndex])
            FaultLenth = length(UpperMatrix[:,1])
            for i=1:FaultLenth; for j=1:FaultLenth; if i >= j; UpperMatrix[i,j] = 0; end; end; end
            push!(Stiffness_U,UpperMatrix)
            Rank_U = [Rank_U;0]
            ElementRange_U = [ElementRange_U; ElementRange_SR[BlockIndex,:]']

            LowerBlockCount += 1
            LowerMatrix = copy(LoadingStiffnessH[BlockIndex])
            FaultLenth = length(LowerMatrix[:,1])
            for i=1:FaultLenth; for j=1:FaultLenth; if i <= j; LowerMatrix[i,j] = 0; end; end; end
            push!(Stiffness_L,LowerMatrix)
            Rank_L = [Rank_L;-1]
            ElementRange_L = [ElementRange_L; ElementRange_SR[BlockIndex,:]']
            push!(BlockSolution,0)
        elseif ElementRange_SR[BlockIndex,1] > ElementRange_SR[BlockIndex,3]
            UpperBlockCount += 1
            push!(Stiffness_U,LoadingStiffnessH[BlockIndex])
            ElementRange_U = [ElementRange_U; ElementRange_SR[BlockIndex,:]']
            Rank_U = [Rank_U;Ranks_Shear[BlockIndex]]
        else
            LowerBlockCount += 1
            push!(Stiffness_L,LoadingStiffnessH[BlockIndex])
            
            ElementRange_L = [ElementRange_L; ElementRange_SR[BlockIndex,:]']
            Rank_L = [Rank_L;Ranks_Shear[BlockIndex]]
            
            if Ranks_Shear[BlockIndex] == 0
                push!(BlockSolution,0)
            else
                push!(BlockSolution,zeros(BlockElementCount))    
            end
        end
    end

    ElementRange_U = ElementRange_U[2:end,:]
    ElementRange_L = ElementRange_L[2:end,:]
    Rank_L = Rank_L[2:end,:]
    Rank_U = Rank_U[2:end,:]

    BlockCount_L = length(ElementRange_L[:,1])
    BlockCount_U = length(ElementRange_U[:,1])
        
    # Lower triangular fault Block-Element relationship
    BlocksHaveFault_NonDiagonal = Any[]
    BlocksHaveFault_Diagonal = Any[]
    for FaultIndex = 1:FaultCount
        BlockHaveThisFault_Diagonal = [0]
        BlockHaveThisFault_NonDiagonal = [0]
        for Blockidx = 1:BlockCount_L
            if ElementRange_L[Blockidx,3] <= FaultIndex && FaultIndex <= ElementRange_L[Blockidx,4] 
                if ElementRange_L[Blockidx,1] == ElementRange_L[Blockidx,3] 
                    BlockHaveThisFault_Diagonal = [BlockHaveThisFault_Diagonal; Blockidx]            
                else
                    BlockHaveThisFault_NonDiagonal= [BlockHaveThisFault_NonDiagonal; Blockidx]     
                end
            end
        end
        BlockHaveThisFault_Diagonal = BlockHaveThisFault_Diagonal[2:end]
        BlockHaveThisFault_NonDiagonal = BlockHaveThisFault_NonDiagonal[2:end]
        push!(BlocksHaveFault_Diagonal,BlockHaveThisFault_Diagonal)
        push!(BlocksHaveFault_NonDiagonal,BlockHaveThisFault_NonDiagonal)
    end



    # calculate Gause-Seidel

    println("Gauss-Seidel")
    Disp_Current = zeros(FaultCount)
    Disp_Prev = zeros(FaultCount)
    MaxDiff = 1
    iteration = 0
    # Par_ElementDivision_U = [0, 2, 5, 40,50,BlockCount_U]
    @time while MaxDiff > DispCritEps
        MaxDiffOld = MaxDiff
        iteration +=1
        if iteration % 20 ==0
        println("Iteration until 1 > ",MaxDiff / DispCritEps)
        end
        
        Disp_PrevStep = copy(Disp_Current)
        Disp_Prev = copy(Disp_Current)
        Disp_Current = zeros(FaultCount)

        Stress_Prev = zeros(FaultCount)
        # Stress_Prev = HmatSolver_Pararllel(Disp_Prev, Stiffness_U, ElementRange_U, Par_ElementDivision_U, ThreadCount, zeros(FaultCount, ThreadCount))
        for Blockidx = 1:BlockCount_U
            
            Stress_Prev[ElementRange_U[Blockidx,3]:ElementRange_U[Blockidx,4]] +=
                Stiffness_U[Blockidx] * Disp_Prev[ElementRange_U[Blockidx,1]:ElementRange_U[Blockidx,2]]    
    
        end

            
        Disp_Current = zeros(FaultCount)
        Disp_Current[1] = (InitialShearStress[1]  +Stress_Prev[1])./ K_Self[1]
        Calculated = zeros(Int,BlockCount_L)
        for FaultIndex = 2:FaultCount
            Stress_Current = 0.0
            for Blockidx in BlocksHaveFault_Diagonal[FaultIndex] # Diagonal Block update
                EleNumInThis = FaultIndex - ElementRange_L[Blockidx,3] + 1 
                Stress_Current += Stiffness_L[Blockidx][EleNumInThis,:]' * Disp_Current[ElementRange_L[Blockidx,1]: ElementRange_L[Blockidx,2]]
                Calculated[Blockidx] = 1
            end    

            for Blockidx in BlocksHaveFault_NonDiagonal[FaultIndex] # Non-diagonal Block update
                EleNumInThis = FaultIndex - ElementRange_L[Blockidx,3] + 1 
                if Calculated[Blockidx] == 0 # calculate mat x vactor if not calculated yet
                    BlockSolution[Blockidx] = Stiffness_L[Blockidx] * Disp_Current[ElementRange_L[Blockidx,1]:ElementRange_L[Blockidx,2]]
                    Calculated[Blockidx] = 1    
                end
                    Stress_Current += BlockSolution[Blockidx][EleNumInThis]
            end 

            StressGap = (InitialShearStress[FaultIndex] + Stress_Current + Stress_Prev[FaultIndex]) 
            Disp_Current[FaultIndex] = (1 - w_factor) * Disp_Prev[FaultIndex] + w_factor * StressGap ./ K_Self[FaultIndex]
            # Disp_Current[FaultIndex] = StressGap ./ K_Self[FaultIndex]
            # Disp_Prev[FaultIndex] = 0.0
        end
            MaxDiff = maximum(abs.(Disp_PrevStep - Disp_Current) )
            # if MaxDiffOld < MaxDiff
            #     w_factor = 1.0
            # end
    end

    println("Iteration GS: ",iteration)

    return Disp_Current
end
