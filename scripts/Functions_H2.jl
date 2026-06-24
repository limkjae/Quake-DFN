
function FindMomDaughter(Block_Range_Level, ClusterCount)

    DaughterElem = zeros(Int,ClusterCount,2)
    DaughterElemVec = []
    DaughterLevelAll = zeros(Int,ClusterCount,2)
    IsSingleMom = zeros(Int,ClusterCount)
    IsLeaf = zeros(Int,ClusterCount)
    LeafClusterIdx = []
    SingleMomClusterIdx = []
    RegularClusterIdx = []
    NonLeafClusterIdx = []
    for ClusterIdx = 1:ClusterCount
        MotherLevel = Block_Range_Level[ClusterIdx,3]
        MotherHierarchy = Block_Range_Level[ClusterIdx,4:3+MotherLevel]
        DaughterCount = 0
        Daughter_i = []
        for DaughterIdx = 1:ClusterIdx
            DaughterLevel = Block_Range_Level[DaughterIdx,3]
            DaughterHierarchy = Block_Range_Level[DaughterIdx,4:3+DaughterLevel-1]
            if MotherLevel == DaughterLevel - 1
                if MotherHierarchy == DaughterHierarchy
                    DaughterCount = DaughterCount + 1
                    DaughterLevelAll[ClusterIdx,DaughterCount] = DaughterLevel
                    DaughterElem[ClusterIdx,DaughterCount] = DaughterIdx
                    push!(Daughter_i, DaughterIdx)
                end
            end
        end
        push!(DaughterElemVec, Daughter_i)
        if DaughterCount == 0
            IsLeaf[ClusterIdx] = 1
            push!(LeafClusterIdx, ClusterIdx)
        elseif DaughterCount == 1 
            IsSingleMom[ClusterIdx] = 1
            push!(SingleMomClusterIdx, ClusterIdx)
            push!(NonLeafClusterIdx, ClusterIdx)
        else 
            push!(RegularClusterIdx, ClusterIdx)
            push!(NonLeafClusterIdx, ClusterIdx)

        end
    end
    return LeafClusterIdx, SingleMomClusterIdx, RegularClusterIdx, NonLeafClusterIdx, IsSingleMom, IsLeaf, DaughterElemVec, DaughterElem
end


function FindW_SN(ClusterCount, Input_Segment, Block_Range_Level, DaughterElemVec, IsLeaf, IsSingleMom)

    Cluster_SkeletonLocal_S = Vector{Int}[]
    Cluster_SkeletonGlobal_S = Vector{Int}[]
    Cluster_W_S =  Matrix{Float64}[]
    
    Cluster_SkeletonLocal_N = Vector{Int}[]
    Cluster_SkeletonGlobal_N = Vector{Int}[]
    Cluster_W_N =  Matrix{Float64}[]

    Cluster_ClusterRange = Vector{Int}[]
    
    SkeletonPrevLevel_S = []
    SkeletonThisLevel_S = []
    SkeletonPrevLevel_N = []
    SkeletonThisLevel_N = []
    SkeletonAll = []
    OrderEachLevel = 1 
    RangeInEachLevel = zeros(Int, ClusterCount,2)
    FaultCenter = Input_Segment[:,1:3]

    for Cluster_i = 1:ClusterCount
        
        ClusterLevel = Block_Range_Level[Cluster_i,3]
        DaughterElem_i = DaughterElemVec[Cluster_i]
        ClusterRange = [Block_Range_Level[Cluster_i,1], Block_Range_Level[Cluster_i,2]]
        println(Cluster_i)
        if Cluster_i > 1 ###### Calculate the order of this cluster in each level
            if ClusterLevel == Block_Range_Level[Cluster_i-1,3] # if current level is same to the previous
                OrderEachLevel = OrderEachLevel + 1
            else ##### if level changes, update the skeletons
                OrderEachLevel = 1
                SkeletonPrevLevel_S = SkeletonThisLevel_S
                SkeletonThisLevel_S = []
                SkeletonPrevLevel_N = SkeletonThisLevel_N
                SkeletonThisLevel_N = []
            end
        end

        RangeInEachLevel[Cluster_i,:] = [(OrderEachLevel-1)*Rank + 1, OrderEachLevel * Rank]

        if IsLeaf[Cluster_i] == 1 # This is Leaf Cluster
            

            ClusterCenter = Block_Ctr_Diam[Cluster_i,1:3] ### Read Ctr and Diam from Full Matrix
            ClusterDiameter = Block_Ctr_Diam[Cluster_i,4]
            ClusterRangeCollect = collect(ClusterRange[1]:ClusterRange[2])

            ProxyElements_i_S, Skeleton_Local_i_S, Skeleton_Global_i_S, W_i_S,
            ProxyElements_i_N, Skeleton_Local_i_N, Skeleton_Global_i_N, W_i_N = 
                SkeletonAndBasis_BuildK_Leaf(Input_Segment, Rank, ClusterCenter, ClusterDiameter, ClusterRangeCollect, 
                                RadMultiple1, RadMultiple2, MinimumProxyElments, RadRange, ShearModulus, PoissonRatio)    
            # println(Cluster_i,"  ", DaughterElem_i, "  ",ClusterRangeCollect)
            SkeletonAll = [SkeletonAll; Skeleton_Global_i_S]

        elseif IsSingleMom[Cluster_i] == 1 # For single mom, just use the Daughter information
            Skeleton_Local_i_S = Cluster_SkeletonLocal_S[DaughterElem_i[1]]
            Skeleton_Global_i_S = Cluster_SkeletonGlobal_S[DaughterElem_i[1]]
            W_i_S = Cluster_W_S[DaughterElem_i[1]]

            Skeleton_Local_i_N = Cluster_SkeletonLocal_N[DaughterElem_i[1]]
            Skeleton_Global_i_N = Cluster_SkeletonGlobal_N[DaughterElem_i[1]]
            W_i_N = Cluster_W_N[DaughterElem_i[1]]

        else #### We only care about the skeletons from the daughter level

            ################# Shear ################
            DaughterSkeleton_S = [Cluster_SkeletonGlobal_S[DaughterElem_i[1]]; Cluster_SkeletonGlobal_S[DaughterElem_i[2]]]
            ClusterCenter_S = vec(sum(FaultCenter[DaughterSkeleton_S,:], dims=1)) / (2*Rank)
            ClusterDiameter_S = maximum(norm.(eachrow(FaultCenter[DaughterSkeleton_S,:] .- ClusterCenter_S'))) *2
            ClusterSkeletonRange_S = [collect((RangeInEachLevel[DaughterElem_i[1],1]) : RangeInEachLevel[DaughterElem_i[1],2])
                                    collect((RangeInEachLevel[DaughterElem_i[2],1]) : RangeInEachLevel[DaughterElem_i[2],2])]
                
            @views ProxyElements_i_S, Skeleton_Local_i_S, Skeleton_LocallyGlobal_i_S, W_i_S = 
            SkeletonAndBasis_BuildK_S(Input_Segment[SkeletonPrevLevel_S,:], Rank, ClusterCenter_S, ClusterDiameter_S, ClusterSkeletonRange_S, 
                                RadMultiple1, RadMultiple2, MinimumProxyElments, RadRange, ShearModulus, PoissonRatio)    

            Skeleton_Global_i_S = DaughterSkeleton_S[Skeleton_Local_i_S]
            ProxyElements_i_S = SkeletonPrevLevel_S[ProxyElements_i_S]

            ################# Normal ################
            DaughterSkeleton_N = [Cluster_SkeletonGlobal_N[DaughterElem_i[1]]; Cluster_SkeletonGlobal_N[DaughterElem_i[2]]]
            ClusterCenter_N = vec(sum(FaultCenter[DaughterSkeleton_N,:], dims=1)) / (2*Rank)
            ClusterDiameter_N = maximum(norm.(eachrow(FaultCenter[DaughterSkeleton_N,:] .- ClusterCenter_N'))) *2
            ClusterSkeletonRange_N = [collect((RangeInEachLevel[DaughterElem_i[1],1]) : RangeInEachLevel[DaughterElem_i[1],2])
                                    collect((RangeInEachLevel[DaughterElem_i[2],1]) : RangeInEachLevel[DaughterElem_i[2],2])]
                
            @views ProxyElements_i_N, Skeleton_Local_i_N, Skeleton_LocallyGlobal_i_N, W_i_N = 
            SkeletonAndBasis_BuildK_N(Input_Segment[SkeletonPrevLevel_N,:], Rank, ClusterCenter_N, ClusterDiameter_N, ClusterSkeletonRange_N, 
                                RadMultiple1, RadMultiple2, MinimumProxyElments, RadRange, ShearModulus, PoissonRatio)    

            Skeleton_Global_i_N = DaughterSkeleton_N[Skeleton_Local_i_N]
            ProxyElements_i_N = SkeletonPrevLevel_N[ProxyElements_i_N]



        end

        SkeletonThisLevel_S = [SkeletonThisLevel_S; Skeleton_Global_i_S]
        SkeletonThisLevel_N = [SkeletonThisLevel_N; Skeleton_Global_i_N]


        push!(Cluster_SkeletonLocal_S, Skeleton_Local_i_S)
        push!(Cluster_SkeletonGlobal_S, Skeleton_Global_i_S)
        push!(Cluster_W_S, W_i_S)   

        push!(Cluster_SkeletonLocal_N, Skeleton_Local_i_N)
        push!(Cluster_SkeletonGlobal_N, Skeleton_Global_i_N)
        push!(Cluster_W_N, W_i_N)   
        push!(Cluster_ClusterRange, ClusterRange)

    end


    return Cluster_SkeletonLocal_S, Cluster_SkeletonGlobal_S, Cluster_W_S, Cluster_ClusterRange,
          Cluster_SkeletonLocal_N, Cluster_SkeletonGlobal_N, Cluster_W_N
end


function FindKd(BlockCount, ElementRange_SR, Admissible, BlockToClusterSR,
                Cluster_SkeletonLocal_S, Cluster_W_S, Cluster_ClusterRange,
                Cluster_SkeletonLocal_N, Cluster_W_N )

    K_d_S = Vector{Matrix{Float64}}(undef, BlockCount)
    K_d_N = Vector{Matrix{Float64}}(undef, BlockCount)

    @time for BlockIdx = 1:BlockCount
        # BlockIdx = 1
        println(BlockIdx)
        BlockRangeS = ElementRange_SR[BlockIdx,1]:ElementRange_SR[BlockIdx,2]
        BlockRangeR = ElementRange_SR[BlockIdx,3]:ElementRange_SR[BlockIdx,4]
        # println(BlockIdx)
        if Admissible[BlockIdx] == 0
            # K_d[BlockIdx], Ki_norm =  StiffnessMatrix_ByParts_Calculation_Rec(Input_Segment[BlockRangeS,:],Input_Segment[BlockRangeR,:],  ShearModulus, PoissonRatio, 1, 1)    
            
            K_d_S[BlockIdx], K_d_N[BlockIdx] =  StiffnessMatrix_ByParts_Calculation_Rec(Input_Segment[BlockRangeS,:],Input_Segment[BlockRangeR,:],  ShearModulus, PoissonRatio, 1, 1)    

            # K_d[BlockIdx] = StiffnessMatrixShear[BlockRangeR,BlockRangeS]
        else
            InteractionPair = [BlockToClusterSR[BlockIdx,:]]
            InteractionPair = convert(Vector{Any},InteractionPair)
            PairIndex = [0]
            PairIdx = 0
            while PairIdx < length(InteractionPair) ## Get the interaction Pair and Pair index
                PairIdx = PairIdx+1
                
                ##### If any of the daughter is 0, calculate just once ######
                DaughterSElements = 0
                DaughterRElements = 0
                DaughterSElements = DaughterElemVec[InteractionPair[PairIdx][1]]
                DaughterRElements = DaughterElemVec[InteractionPair[PairIdx][2]]

    
                ##### Calculate the interaction pair ######
                for DaughterElemS in DaughterSElements
                    for DaughterElemR in DaughterRElements
                        push!(InteractionPair, [DaughterElemS, DaughterElemR])
                        push!(PairIndex, PairIdx)
                    end
                end
            end

            ## Lower level blocks defines 
            LowerLevelBlocks = Vector{Any}(undef, maximum(PairIndex))
            BlockPairCount = zeros(Int, maximum(PairIndex))
            BlockPairCount = convert(Vector{Any},BlockPairCount)
            for Pi in eachindex(PairIndex)
                if PairIndex[Pi] > 0
                    BlockPairCount[PairIndex[Pi]] = BlockPairCount[PairIndex[Pi]]  + 1 
                    if BlockPairCount[PairIndex[Pi]] == 1
                        LowerLevelBlocks[PairIndex[Pi]] = [Pi]
                    else
                        push!(LowerLevelBlocks[PairIndex[Pi]], Pi)
                    end
                end
            end


            ###### InteractionPair  --> information of blocks
            ###### LowerLevelBlocks --> information of daughter block set

            TotalInteractionPair = length(InteractionPair)
            K_id_S = Vector{Matrix{Float64}}(undef, TotalInteractionPair)
            K_id_N = Vector{Matrix{Float64}}(undef, TotalInteractionPair)
            for ClusterPair_i = TotalInteractionPair:-1:1
                # println(ClusterPair_i)
                SIdx = InteractionPair[ClusterPair_i][1]
                RIdx = InteractionPair[ClusterPair_i][2]

                if isassigned(LowerLevelBlocks, ClusterPair_i)

                    if length(LowerLevelBlocks[ClusterPair_i]) == 4
                        # println(ClusterPair_i, "here 1")
                        PairIdx = LowerLevelBlocks[ClusterPair_i]
                        # println(ClusterPair_i," ",PairIdx," ", IsLeaf[SIdx]," ", IsLeaf[RIdx])
                        # println(IsLeaf(SIdx)," ", IsLeaf(RIdx))
                        K_i_S  = [K_id_S[PairIdx[1]] K_id_S[PairIdx[3]]; K_id_S[PairIdx[2]] K_id_S[PairIdx[4]]]                    
                        K_id_S[ClusterPair_i] =  InterpMatrix(K_i_S, Cluster_W_S[SIdx],  Cluster_W_S[RIdx], Cluster_SkeletonLocal_S[SIdx], 
                                    Cluster_SkeletonLocal_S[RIdx]) 

                        K_i_N  = [K_id_N[PairIdx[1]] K_id_N[PairIdx[3]]; K_id_N[PairIdx[2]] K_id_N[PairIdx[4]]]
                        K_id_N[ClusterPair_i] =  InterpMatrix_N(K_i_N, Cluster_W_N[SIdx],  Cluster_W_N[RIdx], Cluster_SkeletonLocal_N[SIdx], 
                                    Cluster_SkeletonLocal_N[RIdx]) 


                    elseif IsSingleMom[SIdx] == 1 && IsSingleMom[RIdx] == 1
                        PairIdx = LowerLevelBlocks[ClusterPair_i]
                            # println("Both SingleMom ",SIdx,"   ",RIdx)
                        
                        K_id_S[ClusterPair_i] = K_id_S[PairIdx[1]]
                        K_id_N[ClusterPair_i] = K_id_N[PairIdx[1]]

                    elseif IsSingleMom[SIdx] == 1
                        PairIdx = LowerLevelBlocks[ClusterPair_i]
                            # println("Source SingleMom ",SIdx)
                        K_i_S  = [K_id_S[PairIdx[1]]; K_id_S[PairIdx[2]]]
                        K_id_S[ClusterPair_i] =  InterpMatrixReceiver(K_i_S,  Cluster_W_S[RIdx], Cluster_SkeletonLocal_S[RIdx])  

                        K_i_N  = [K_id_N[PairIdx[1]]; K_id_N[PairIdx[2]]]
                        K_id_N[ClusterPair_i] =  InterpMatrixReceiver_N(K_i_N,  Cluster_W_N[RIdx], Cluster_SkeletonLocal_N[RIdx])  

                        # println(ClusterPair_i, "here 1")
                    elseif IsSingleMom[RIdx] == 1
                            # println("Receiver SingleMom",RIdx)
                        PairIdx = LowerLevelBlocks[ClusterPair_i]
                        K_i_S  = [K_id_S[PairIdx[1]] K_id_S[PairIdx[2]]]
                        K_id_S[ClusterPair_i] =  InterpMatrixSource(K_i_S,  Cluster_W_S[SIdx], Cluster_SkeletonLocal_S[SIdx])  

                        K_i_N  = [K_id_N[PairIdx[1]] K_id_N[PairIdx[2]]]
                        K_id_N[ClusterPair_i] =  InterpMatrixSource_N(K_i_N,  Cluster_W_N[SIdx], Cluster_SkeletonLocal_N[SIdx])  
                        # println(ClusterPair_i, "here 1")
                    else
                        println(ClusterPair_i, "  ", length(LowerLevelBlocks[ClusterPair_i]))
                        println("Error Here!")
                    end
                
                else ###### This is Leaf-Leaf level ######
                    # println("Leaf-leaf")
                    RangeR = Cluster_ClusterRange[RIdx][1]:Cluster_ClusterRange[RIdx][2]
                    RangeS =  Cluster_ClusterRange[SIdx][1]:Cluster_ClusterRange[SIdx][2]
                    K_i_S, K_i_N =  StiffnessMatrix_ByParts_Calculation_Rec(Input_Segment[RangeS,:],Input_Segment[RangeR,:],  ShearModulus, PoissonRatio, 1, 1)    

                    K_id_S[ClusterPair_i] =  InterpMatrix(K_i_S, Cluster_W_S[SIdx],  Cluster_W_S[RIdx], Cluster_SkeletonLocal_S[SIdx], 
                                                            Cluster_SkeletonLocal_S[RIdx]) 
                    K_id_N[ClusterPair_i] =  InterpMatrix_N(K_i_N, Cluster_W_N[SIdx],  Cluster_W_N[RIdx], Cluster_SkeletonLocal_N[SIdx], 
                                                        Cluster_SkeletonLocal_N[RIdx]) 
                        # println(ClusterPair_i, "here 5")
                end
            end
            K_d_S[BlockIdx] = K_id_S[1]
            K_d_N[BlockIdx] = K_id_N[1]
        end


        # if BlockIdx > 500
            if BlockIdx % 500 == 0
               print("\033c")   
            end
        # end
    end

return K_d_S, K_d_N
end






function SkeletonAndBasis2(K, Rank, ClusterCenter, ClusterDiam, ClusterRange, FaultCenter, 
    RadMultiple1, RadMultiple2, MinimumProxyElments, RadRange)

    ProxyElements = []
    FaultCount_i  = length(FaultCenter[:,1])
    while length(ProxyElements) < MinimumProxyElments
        ProxyElements = []
        for i=1:FaultCount_i
            FaultDistance = norm(FaultCenter[i,:] - ClusterCenter)
            if ClusterDiam * RadMultiple1 - RadRange <  FaultDistance < ClusterDiam * RadMultiple1 + RadRange ||
                ClusterDiam * RadMultiple2 - RadRange <  FaultDistance < ClusterDiam * RadMultiple2 + RadRange
                ProxyElements = [ProxyElements; i]
            end
        end
        RadRange = RadRange + 100
        
        # RadMultiple2 = RadMultiple2 * 0.7
    end
    @views SK_Proxy_Matrix = K[ClusterRange,ProxyElements]
    Skeleton_Local = idfact(SK_Proxy_Matrix', rank = Rank)[:sk]
    # println(length(Skeleton_Local))
    Skeleton_Global =  Skeleton_Local .+ ClusterRange[1] .- 1
    W = SK_Proxy_Matrix'[:,Skeleton_Local] \ SK_Proxy_Matrix'

    return ProxyElements, Skeleton_Local, Skeleton_Global, W
end

function SkeletonAndBasis_BuildK_Leaf(Input_Segment, Rank, ClusterCenter, ClusterDiam, ClusterRange,  
    RadMultiple1, RadMultiple2, MinimumProxyElments, RadRange, ShearModulus, PoissonRatio)   

    FaultCenter = Input_Segment[:,1:3]
    ProxyElements = []
    FaultCount_i  = length(FaultCenter[:,1])
    while length(ProxyElements) < MinimumProxyElments
        ProxyElements = []
        for i=1:FaultCount_i
            FaultDistance = norm(FaultCenter[i,:] - ClusterCenter)
            if ClusterDiam * RadMultiple1 - RadRange <  FaultDistance < ClusterDiam * RadMultiple1 + RadRange ||
                ClusterDiam * RadMultiple2 - RadRange <  FaultDistance < ClusterDiam * RadMultiple2 + RadRange
                ProxyElements = [ProxyElements; i]
            end
        end
        RadRange = RadRange + 100
        
        # RadMultiple2 = RadMultiple2 * 0.7
    end
    
    SK_Proxy_Matrix_S, SK_Proxy_Matrix_N = 
    StiffnessMatrix_ByParts_Calculation_Rec(Input_Segment[ProxyElements,:],Input_Segment[ClusterRange,:],  ShearModulus, PoissonRatio,
                                        1, 1)    

    Skeleton_Local_S = idfact(SK_Proxy_Matrix_S', rank = Rank)[:sk]
    Skeleton_Global_S =  Skeleton_Local_S .+ ClusterRange[1] .- 1
    W_S = SK_Proxy_Matrix_S'[:,Skeleton_Local_S] \ SK_Proxy_Matrix_S'



    Skeleton_Local_N = idfact(SK_Proxy_Matrix_N', rank = Rank)[:sk]
    RankRange = size(SK_Proxy_Matrix_N)[1]    
    if length(Skeleton_Local_N) != Rank # if Skeleton_Local is smaller than Rank * 2, add some random number.
        # println(SK_Proxy_Matrix)
        while length(Skeleton_Local_N) < Rank
            # println("Adding Random Skeleton")
            RandElem = rand(1:RankRange)
            if RandElem in Skeleton_Local_N
            else
                push!(Skeleton_Local_N, RandElem)
            end
        end        
    end

    Skeleton_Global_N =  Skeleton_Local_N .+ ClusterRange[1] .- 1
    W_N = SK_Proxy_Matrix_N'[:,Skeleton_Local_N] \ SK_Proxy_Matrix_N' 

    ProxyElements_S = copy(ProxyElements)
    ProxyElements_N =  copy(ProxyElements)

    return ProxyElements_S, Skeleton_Local_S, Skeleton_Global_S, W_S,
           ProxyElements_N, Skeleton_Local_N, Skeleton_Global_N, W_N
end

function SkeletonAndBasis_BuildK_S(Input_Segment, Rank, ClusterCenter, ClusterDiam, ClusterRange,  
    RadMultiple1, RadMultiple2, MinimumProxyElments, RadRange, ShearModulus, PoissonRatio)   

    FaultCenter = Input_Segment[:,1:3]
    ProxyElements = []
    FaultCount_i  = length(FaultCenter[:,1])
    while length(ProxyElements) < MinimumProxyElments
        ProxyElements = []
        for i=1:FaultCount_i
            FaultDistance = norm(FaultCenter[i,:] - ClusterCenter)
            if ClusterDiam * RadMultiple1 - RadRange <  FaultDistance < ClusterDiam * RadMultiple1 + RadRange ||
                ClusterDiam * RadMultiple2 - RadRange <  FaultDistance < ClusterDiam * RadMultiple2 + RadRange
                ProxyElements = [ProxyElements; i]
            end
        end
        RadRange = RadRange + 100
        
        # RadMultiple2 = RadMultiple2 * 0.7
    end
    
    SK_Proxy_Matrix, SK_Proxy_MatrixNorm = 
    StiffnessMatrix_ByParts_Calculation_Rec(Input_Segment[ProxyElements,:],Input_Segment[ClusterRange,:],  ShearModulus, PoissonRatio,
                                        1, 1)    
# println(ClusterRange[1])
    # @views SK_Proxy_Matrix = K[ClusterRange,ProxyElements]
    Skeleton_Local = idfact(SK_Proxy_Matrix', rank = Rank)[:sk]
    # println(length(Skeleton_Local))
    Skeleton_Global =  Skeleton_Local .+ ClusterRange[1] .- 1
    W = SK_Proxy_Matrix'[:,Skeleton_Local] \ SK_Proxy_Matrix'

    return ProxyElements, Skeleton_Local, Skeleton_Global, W
end

function SkeletonAndBasis_BuildK_N(Input_Segment, Rank, ClusterCenter, ClusterDiam, ClusterRange,  
    RadMultiple1, RadMultiple2, MinimumProxyElments, RadRange, ShearModulus, PoissonRatio)   

    FaultCenter = Input_Segment[:,1:3]
    ProxyElements = []
    FaultCount_i  = length(FaultCenter[:,1])
    while length(ProxyElements) < MinimumProxyElments
        ProxyElements = []
        for i=1:FaultCount_i
            FaultDistance = norm(FaultCenter[i,:] - ClusterCenter)
            if ClusterDiam * RadMultiple1 - RadRange <  FaultDistance < ClusterDiam * RadMultiple1 + RadRange ||
                ClusterDiam * RadMultiple2 - RadRange <  FaultDistance < ClusterDiam * RadMultiple2 + RadRange
                ProxyElements = [ProxyElements; i]
            end
        end
        RadRange = RadRange + 100
        
        # RadMultiple2 = RadMultiple2 * 0.7
    end
    
    SK_Proxy_MatrixShear, SK_Proxy_Matrix = 
    StiffnessMatrix_ByParts_Calculation_Rec(Input_Segment[ProxyElements,:],Input_Segment[ClusterRange,:],  ShearModulus, PoissonRatio,
                                        1, 1)  

    Skeleton_Local = idfact(SK_Proxy_Matrix', rank = Rank)[:sk]

    RankRange = size(SK_Proxy_Matrix)[1]    
    if length(Skeleton_Local) != Rank # if Skeleton_Local is smaller than Rank * 2, add some random number.
        # println(SK_Proxy_Matrix)
        println("Adding Random Skeleton")
        while length(Skeleton_Local) < Rank
            RandElem = rand(1:RankRange)
            if RandElem in Skeleton_Local
            else
                push!(Skeleton_Local, RandElem)
            end
        end        
    end
    

    Skeleton_Global =  Skeleton_Local .+ ClusterRange[1] .- 1
    W = SK_Proxy_Matrix'[:,Skeleton_Local] \ SK_Proxy_Matrix'

    return ProxyElements, Skeleton_Local, Skeleton_Global, W

end



function SkeletonAndBasis_BuildK_Normal(Input_Segment, Rank, ClusterCenter, ClusterDiam, ClusterRange,  
    RadMultiple1, RadMultiple2, MinimumProxyElments, RadRange, ShearModulus, PoissonRatio)   

    FaultCenter = Input_Segment[:,1:3]
    ProxyElements = []
    FaultCount_i  = length(FaultCenter[:,1])
    while length(ProxyElements) < MinimumProxyElments
        ProxyElements = []
        for i=1:FaultCount_i
            FaultDistance = norm(FaultCenter[i,:] - ClusterCenter)
            if ClusterDiam * RadMultiple1 - RadRange <  FaultDistance < ClusterDiam * RadMultiple1 + RadRange ||
                ClusterDiam * RadMultiple2 - RadRange <  FaultDistance < ClusterDiam * RadMultiple2 + RadRange
                ProxyElements = [ProxyElements; i]
            end
        end
        RadRange = RadRange + 100
        
        # RadMultiple2 = RadMultiple2 * 0.7
    end
    
    SK_Proxy_MatrixShear, SK_Proxy_Matrix = 
    StiffnessMatrix_ByParts_Calculation_Rec(Input_Segment[ProxyElements,:],Input_Segment[ClusterRange,:],  ShearModulus, PoissonRatio,
                                        1, 1)  

    Skeleton_Local = idfact(SK_Proxy_Matrix', rank = Rank)[:sk]
    println(maximum(SK_Proxy_Matrix) - minimum(SK_Proxy_Matrix), "   ", length(Skeleton_Local))
    
    println(size(SK_Proxy_Matrix),  length(Skeleton_Local), "  ",ProxyElements)

    RankRange = size(SK_Proxy_Matrix)[1]    
    if length(Skeleton_Local) != Rank # if Skeleton_Local is smaller than Rank * 2, add some random number.
        # println(SK_Proxy_Matrix)
        while length(Skeleton_Local) < Rank
            RandElem = rand(1:RankRange)
            if RandElem in Skeleton_Local
            else
                push!(Skeleton_Local, RandElem)
            end
        end        
    end
    println(size(SK_Proxy_Matrix),  length(Skeleton_Local), "  ",ProxyElements)

    Skeleton_Global =  Skeleton_Local .+ ClusterRange[1] .- 1
    W = SK_Proxy_Matrix'[:,Skeleton_Local] \ SK_Proxy_Matrix'

    return ProxyElements, Skeleton_Local, Skeleton_Global, W
end

function InterpMatrix(K, W_s,W_r, SkS, SkR)
    # Interpolate the matrix to the basis of each cluster for each source and receiver
    # K: Original Matrix, W: W value (basis), Sk: Skeleton Sorce and Receiver
    
    W_i = W_s
    U_i = W_r

    Interp_S = K[:,SkS] \ K
    T_s = Interp_S / W_i

    Interp_R = K'[:,SkR] \ K'
    T_r = Interp_R / U_i
    K_id = T_r' * K[SkR,SkS] * T_s
    
    return K_id
end

function InterpMatrix_N(K, W_s,W_r, SkS, SkR)
    # Interpolate the matrix to the basis of each cluster for each source and receiver
    # K: Original Matrix, W: W value (basis), Sk: Skeleton Sorce and Receiver
    
    W_i = W_s
    U_i = W_r

    Interp_S = pinv(K[:, SkS]) * K
    T_s = Interp_S / W_i

    Interp_R = pinv(K'[:, SkR]) * K'
    T_r = Interp_R / U_i
    K_id = T_r' * K[SkR,SkS] * T_s


    return K_id
end



function InterpMatrixReceiver(K, W_r, SkR)
    # Interpolate the matrix to the basis of each cluster for each source and receiver
    # K: Original Matrix, W: W value (basis), Sk: Skeleton Sorce and Receiver
    
    U_i = W_r

    Interp_R = K'[:,SkR] \ K'
    T_r = Interp_R / U_i
    K_id = T_r' * K[SkR,:] 

    return K_id
end
    
function InterpMatrixReceiver_N(K, W_r, SkR)
    # Interpolate the matrix to the basis of each cluster for each source and receiver
    # K: Original Matrix, W: W value (basis), Sk: Skeleton Sorce and Receiver
    
    U_i = W_r

    # Interp_R = K'[:,SkR] \ K'
    # Interp_R = qr(K'[:,SkR] , Val(true)) \ K'
    Interp_R = pinv(K'[:, SkR]) * K'
    T_r = Interp_R / U_i
    K_id = T_r' * K[SkR,:] 



    return K_id
end

function InterpMatrixSource(K, W_s, SkS)
    # Interpolate the matrix to the basis of each cluster for each source and receiver
    # K: Original Matrix, W: W value (basis), Sk: Skeleton Sorce and Receiver
    
    W_i = W_s
    
    # println(inv(K))

    Interp_S = K[:,SkS] \ K
    # Interp_S = pinv(K[:, SkS]) * K
    # Interp_S = qr(K[:, SkS], Val(true)) \ K
    T_s = Interp_S / W_i
    K_id = K[:,SkS] * T_s

    return K_id
end

function InterpMatrixSource_N(K, W_s, SkS)
    # Interpolate the matrix to the basis of each cluster for each source and receiver
    # K: Original Matrix, W: W value (basis), Sk: Skeleton Sorce and Receiver
    
    W_i = W_s
    
    # println(inv(K))
    # if cond(K[:, SkS]) > 1e10
    # else 
        # Interp_S = K[:,SkS] \ K
    # end
    # Interp_S = K[:,SkS] \ K
    # Interp_S = pinv(K[:, SkS]) * K
    # Interp_S = qr(K[:, SkS], Val(true)) \ K
    Interp_S = pinv(K[:, SkS]) * K
    T_s = Interp_S / W_i
    K_id = K[:,SkS] * T_s

    return K_id
end





function MatVecMulti_Parallel_Inadmissible!(Stress::Vector{Float64}, InAdmissibleBlockSep::Vector{Vector{Int64}}, ElementRange_SR::Matrix{Int64},
                             K_d::Vector{Matrix{Float64}}, Disp::Vector{Float64}, ThreadCount::Int64,  
                             thread_buffers::Vector{Vector{Float64}})::Vector{Float64}
    
    # @views Threads.@threads 
    # @sync begin
    #     @views Threads.@threads  for Ti = 1:ThreadCount
    #         for BlockIdx in InAdmissibleBlockSep[Ti]
    #             BlockRangeS = ElementRange_SR[BlockIdx,1]:ElementRange_SR[BlockIdx,2]
    #             BlockRangeR = ElementRange_SR[BlockIdx,3]:ElementRange_SR[BlockIdx,4]
    #             Stress[BlockRangeR] += K_d[BlockIdx] * Disp[BlockRangeS]

    #         end
    #     end
    # end

    @sync Threads.@threads for Ti = 1:ThreadCount
        local_stress = thread_buffers[Ti]
        
        for BlockIdx in InAdmissibleBlockSep[Ti]
            # Create views to avoid copying memory
            SR_row = @view ElementRange_SR[BlockIdx, :]
            BlockRangeS = SR_row[1]:SR_row[2]
            BlockRangeR = SR_row[3]:SR_row[4]
            
            # Views of the specific chunks we are working on
            stress_view = @view local_stress[BlockRangeR]
            disp_view   = @view Disp[BlockRangeS]
            K_matrix    = K_d[BlockIdx]
            
            mul!(stress_view, K_matrix, disp_view, 1.0, 1.0)
        end
    end

    for Ti = 1:ThreadCount
        Stress .+= thread_buffers[Ti]
    end

    return Stress
end

function MatVecMulti_Parallel_InadmissibleP2(Stress::Vector{Float64}, InAdmissibleBlockSep::Vector{Vector{Int64}}, ElementRange_SR::Matrix{Int64},
                             K_d::Vector{Matrix{Float64}}, Disp::Vector{Float64}, ThreadCount::Int64,  
                             thread_buffers::Vector{Vector{Float64}}, Stress_EachThread::Vector{Any}, 
                             Stress_RangeEachThread::Matrix{Int64})::Vector{Float64}
    

# Zero out Stress before accumulation if this functions handles the full initialization
    fill!(Stress, 0.0)

    Threads.@threads for Ti = 1:ThreadCount
        blocks = InAdmissibleBlockSep[Ti]
        
        for i in eachindex(blocks)
            BlockIdx = blocks[i]
            
            # Extract bounds directly to avoid row-view overhead
            @inbounds r_start = ElementRange_SR[BlockIdx, 1]
            @inbounds r_end   = ElementRange_SR[BlockIdx, 2]
            @inbounds r_three = ElementRange_SR[BlockIdx, 3]
            @inbounds r_four  = ElementRange_SR[BlockIdx, 4]
            
            BlockRangeS = r_start:r_end
            # If BlockRangeR maps directly to global indices in Stress:
            BlockRangeR = r_three:r_four 
            
            # Subviews directly into global arrays
            stress_view = @view Stress[BlockRangeR]
            disp_view   = @view Disp[BlockRangeS]
            K_matrix    = K_d[BlockIdx]
            
            # Direct in-place accumulation safely isolated by thread
            mul!(stress_view, K_matrix, disp_view, 1.0, 1.0)
        end
    end

    return Stress
end



function InteractionBlocks!(AdmissibleBlockSep::Vector{Vector{Int64}}, BlockToClusterSR::Matrix{Int64},
                            Sig_Cluster::Vector{Vector{Float64}}, D_Cluster::Vector{Vector{Float64}}, 
                            K_d::Vector{Matrix{Float64}}, ThreadCount::Int64)

    cluster_locks = [ReentrantLock() for _ in eachindex(Sig_Cluster)]

    @sync Threads.@threads for Ti = 1:ThreadCount
         for BlockIdx in AdmissibleBlockSep[Ti] 
            r_idx = BlockToClusterSR[BlockIdx, 2]
            s_idx = BlockToClusterSR[BlockIdx, 1]
            
            sig_view  = Sig_Cluster[r_idx]
            disp_view = D_Cluster[s_idx]
            K_matrix  = K_d[BlockIdx]
            
            # Safely lock only the cluster being updated
            lock(cluster_locks[r_idx]) do
                mul!(sig_view, K_matrix, disp_view, 1.0, 1.0)
            end
        end
    end
    return Sig_Cluster
end

function ForwardCalculation(ClusterCount::Int64, IsLeaf::Vector{Int64}, IsSingleMom::Vector{Int64},
                        Cluster_ClusterRange::Vector{Vector{Int64}}, D_Cluster::Vector{Vector{Float64}},
                        Cluster_W::Vector{Matrix{Float64}}, DaughterElem::Matrix{Int64}, Disp::Vector{Float64})

    @inbounds for ClusterIdx = 1:ClusterCount
        
        if IsLeaf[ClusterIdx] == 1
            r_start = Cluster_ClusterRange[ClusterIdx][1]
            r_end   = Cluster_ClusterRange[ClusterIdx][2]
            disp_view = @view Disp[r_start:r_end]            
            mul!(D_Cluster[ClusterIdx], Cluster_W[ClusterIdx], disp_view)
            
        elseif IsSingleMom[ClusterIdx] == 1
            DaughterCluster = DaughterElem[ClusterIdx, 1]
            D_Cluster[ClusterIdx] = D_Cluster[DaughterCluster]
            
        else
            combined_view =  [D_Cluster[DaughterElem[ClusterIdx,1]]; D_Cluster[DaughterElem[ClusterIdx,2]]]
            mul!(D_Cluster[ClusterIdx], Cluster_W[ClusterIdx], combined_view)
        end
    end
    return D_Cluster

end


function BackwardCalculation(NonLeafClusterIndex::Vector{Any}, LeafClusterIndex::Vector{Any},
                            IsSingleMom::Vector{Int64}, DaughterElem::Matrix{Int64}, 
                            Sig_Cluster::Vector{Vector{Float64}}, Rank::Int64,
                            Block_Range_Level::Matrix{Int64}, Stress::Vector{Float64},
                            Cluster_W::Vector{Matrix{Float64}},ThreadCount::Int64)


    # --- 1. PRE-ALLOCATION FOR MULTI-THREADING ---
    # To use mul! without race conditions or allocations, give each thread its own buffer.
    # If Rank is small (e.g., <= 32), consider StaticArrays instead of this.
    ntreads = ThreadCount
    buffers = [Vector{Float64}(undef, 2 * Rank) for _ in 1:ntreads]

    # --- 2. DOWNWARD PASS (NON-LEAF) ---
    # WARNING: Cannot easily parallelize this loop with @threads because parent clusters 
    # push data down to daughters. If a level has dependencies, threading will cause data races.
    # We use @inbounds and mul! to maximize sequential performance.
    @inbounds for ClusterIdx in reverse(NonLeafClusterIndex)
        if IsSingleMom[ClusterIdx] == 1
            Dtr1 = DaughterElem[ClusterIdx, 1]
            Sig_Cluster[Dtr1] .+= Sig_Cluster[ClusterIdx]
        else
            Dtr1 = DaughterElem[ClusterIdx, 1]
            Dtr2 = DaughterElem[ClusterIdx, 2]
            
            # Non-allocating matrix-vector multiplication using thread-safe buffer
            tid = threadid()
            buf = buffers[tid]
            
            # buf = Cluster_W' * Sig_Cluster
            mul!(buf, Cluster_W[ClusterIdx]', Sig_Cluster[ClusterIdx])
            
            # View-based in-place accumulation
            Sig_Cluster[Dtr1] .+= @view buf[1:Rank]
            Sig_Cluster[Dtr2] .+= @view buf[Rank+1:2*Rank]
        end
    end

    # --- 3. LEAF PASS (PARALLELIZED) ---
    # Safe to parallelize if different Leaf clusters write to independent slices of `Stress`.
    # If leaf ranges overlap, change `@threads` to a sequential loop or use locks.
    @threads for ClusterIdx in reverse(LeafClusterIndex)
        @inbounds begin
            r_start = Block_Range_Level[ClusterIdx, 1]
            r_end   = Block_Range_Level[ClusterIdx, 2]
            
            # Extract views to avoid allocation
            stress_slice = @view Stress[r_start:r_end]
            W_T = Cluster_W[ClusterIdx]'
            Sig = Sig_Cluster[ClusterIdx]
            
            # In-place fused multiply-add: stress_slice += W_T * Sig
            mul!(stress_slice, W_T, Sig, 1.0, 1.0)
        end
    end



    return Stress
end


function H2MatVecCalculation(K_d, Cluster_W, Stress, InAdmissibleBlockSep, ElementRange_SR,  
                            Disp, ThreadCount, thread_buffers,
                            ClusterCount, IsLeaf, IsSingleMom, Cluster_ClusterRange, D_Cluster, 
                            DaughterElem, AdmissibleBlockSep, BlockToClusterSR,  Sig_Cluster, 
                            NonLeafClusterIndex, LeafClusterIndex, Rank, Block_Range_Level)

    MatVecMulti_Parallel_Inadmissible!(Stress, InAdmissibleBlockSep, ElementRange_SR, K_d, Disp, ThreadCount, thread_buffers)     
    
    # @time Stress = MatVecMulti_Parallel_InadmissibleP2(Stress, InAdmissibleBlockSep, ElementRange_SR, K_d, Disp,
    #                                 ThreadCount, thread_buffers, Stress_EachThread, Stress_RangeEachThread)    
    

    ForwardCalculation(ClusterCount, IsLeaf, IsSingleMom, Cluster_ClusterRange, D_Cluster,
                            Cluster_W, DaughterElem, Disp)


    InteractionBlocks!(AdmissibleBlockSep, BlockToClusterSR,  Sig_Cluster, D_Cluster, 
                                K_d,ThreadCount)


    BackwardCalculation(NonLeafClusterIndex, LeafClusterIndex, IsSingleMom, DaughterElem, 
                                Sig_Cluster, Rank, Block_Range_Level, Stress, Cluster_W, ThreadCount) 


    return Stress
end


function Ad_InAd_Blocks_EachThreads(BlockCount, Admissible, ThreadCount, FaultCount, 
                                Block_Range_Level, BlockToClusterSR)

    AdmissibleBlock = Int[]
    InAdmissibleBlock = Int[]

    for BlockIdx = 1:BlockCount
        if Admissible[BlockIdx] == 1 
            push!(AdmissibleBlock, BlockIdx)
        else
            push!(InAdmissibleBlock, BlockIdx)
        end
    end

    InAdmissibleBlockSep = Vector{Vector{Int64}}(undef, ThreadCount)
        SepUnit = ceil(Int, FaultCount / ThreadCount)
    for i=1:ThreadCount; InAdmissibleBlockSep[i] = []; end
        for BlockIndex in InAdmissibleBlock
            WhereToPutIdx = ceil(Int, Block_Range_Level[BlockToClusterSR[BlockIndex,2],2]/SepUnit)
            push!(InAdmissibleBlockSep[WhereToPutIdx], BlockIndex)
    end

    AdmissibleBlockSep = Vector{Vector{Int64}}(undef, ThreadCount)
    for i=1:ThreadCount; AdmissibleBlockSep[i] = []; end
        AdBlockCount = length(AdmissibleBlock)
        SepUnit = ceil(Int, AdBlockCount/ ThreadCount)
        for BlockIndex = 1:AdBlockCount
            WhereToPutIdx = ceil(Int, BlockIndex/SepUnit)
            push!(AdmissibleBlockSep[WhereToPutIdx], AdmissibleBlock[BlockIndex])
    end

    return AdmissibleBlockSep, InAdmissibleBlockSep
end
