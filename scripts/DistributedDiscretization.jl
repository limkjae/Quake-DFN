

using DelimitedFiles
# using Base
using PyPlot
using PyCall
using JLD2
using LowRankApprox
using Clustering
using LinearAlgebra
# using Distributed
@pyimport matplotlib.patches as patches

pygui(true)

HowManyDistribution = 5

include("../scripts/Functions_BuildInputFile.jl")
include("../scripts/Functions_OKADA3D.jl")
include("../scripts/Functions_Plot.jl")
include("../scripts/Functions_Hmatrix.jl")
include("../scripts/Functions_TDstressHS.jl")

isdir("scripts/temp_Discretization") || mkdir("scripts/temp_Discretization")


InputBulkFileName="Input_BulkFaultGeometry.txt"
HMatrixStructureFile = "Input_HmatrixStructure.jld2"

# InputBulkFileName="../Input_BulkFaultGeometry.txt"
# OutputFileName="Input_Discretized.jld2"
# HMatrixStructureFile = "../Input_HmatrixStructure.jld2"

Block_Ctr_Diam, Block_Range_Level, Input_Segment, LoadingFaultExist, LoadingFaultCount,
MinimumElementsToCut, ArrangePoint, Admissible, ElementRange_SR, Switch_StrikeSlip_or_ReverseNormal, 
ShearModulus, PoissonRatio, RockDensity, DropCrit, DropCritNormalStressMultiplier, MinimumNS, RorT =
    load(HMatrixStructureFile, "Block_Ctr_Diam", "Block_Range_Level", "Input_Segment", "LoadingFaultExist", "LoadingFaultCount",
        "MinimumElementsToCut", "ArrangePoint", "Admissible", "ElementRange_SR", "Switch_StrikeSlip_or_ReverseNormal", 
        "ShearModulus", "PoissonRatio", "RockDensity", "DropCrit", "DropCritNormalStressMultiplier", "MinimumNS", "RorT")


BlockElementCount = (ElementRange_SR[:,2] - ElementRange_SR[:,1] .+ 1) .* (ElementRange_SR[:,4] - ElementRange_SR[:,3] .+ 1)  

BlockElementIncrementNorm = cumsum(BlockElementCount) / sum(BlockElementCount)
DivisionCut = zeros(Int, HowManyDistribution)
DivisionBlockItoF = zeros(Int, HowManyDistribution,2)
TotalBlock = length(BlockElementIncrementNorm)


for BlockIndex in eachindex(BlockElementIncrementNorm)
    for DivisionIndex = 1 : HowManyDistribution - 1
        DivisionCrit =  DivisionIndex / HowManyDistribution 
        if BlockElementIncrementNorm[BlockIndex] < DivisionCrit
            DivisionCut[DivisionIndex] = BlockIndex
           
        end
    end
end
DivisionCut[end] = length(BlockElementIncrementNorm)

for index = 1:HowManyDistribution
    if index == 1
        DivisionBlockItoF[1,:] = [1, DivisionCut[index]] 
    else
        DivisionBlockItoF[index,:] = [DivisionCut[index-1]+1, DivisionCut[index]] 

    end
end

for DistributeIndex = 1:HowManyDistribution

CodeScript = """

    using DelimitedFiles
    using Base
    using PyPlot
    using PyCall
    using JLD2
    using LowRankApprox
    using Clustering
    using LinearAlgebra

    include("../Functions_BuildInputFile.jl")
    include("../Functions_OKADA3D.jl")
    include("../Functions_Plot.jl")
    include("../Functions_Hmatrix.jl")
    




    function Discritize()

        BlockI = $(DivisionBlockItoF[DistributeIndex,1])
        BlockF = $(DivisionBlockItoF[DistributeIndex,2])

        InputBulkFileName="Input_BulkFaultGeometry.txt"
        OutputFileName="scripts/temp_Discretization/HMatPart_$(DistributeIndex).jld2"
        HMatrixStructureFile = "Input_HmatrixStructure.jld2"

        ##########################       Hmatrix       #############################
        Tolerance = 1e3 #  Hmatrix compression Tolerance. pascal for 1m slip (More approximaion for higher Tolerance). 

        println("H Matrix will be built")
        Block_Ctr_Diam, Block_Range_Level, Input_Segment, LoadingFaultExist, LoadingFaultCount,
        MinimumElementsToCut, ArrangePoint, Admissible, ElementRange_SR, Switch_StrikeSlip_or_ReverseNormal, 
        ShearModulus, PoissonRatio, RockDensity, DropCrit, DropCritNormalStressMultiplier, MinimumNS, RorT =
        load(HMatrixStructureFile, "Block_Ctr_Diam", "Block_Range_Level", "Input_Segment", "LoadingFaultExist", "LoadingFaultCount",
            "MinimumElementsToCut", "ArrangePoint", "Admissible", "ElementRange_SR", "Switch_StrikeSlip_or_ReverseNormal", 
            "ShearModulus", "PoissonRatio", "RockDensity", "DropCrit", "DropCritNormalStressMultiplier", "MinimumNS", "RorT")

        Admissible = Admissible[BlockI:BlockF]
        ElementRange_SR = ElementRange_SR[BlockI:BlockF,:]



        FaultCount = length(Input_Segment[:,1])


        ######################## Read Segmented Input File #######################
        FaultCenter, Fault_a, Fault_b, Fault_Dc, Fault_Theta_i, Fault_V_i, Fault_Friction_i, Fault_NormalStress, 
        Fault_V_Const, Fault_BulkIndex, FaultLengthStrike, FaultLengthDip, FaultStrikeAngle, 
        FaultDipAngle, FaultRakeAngle, FaultLengthStrike_Bulk, FaultLengthDip_Bulk, NormalStiffnessZero = 
            ReadSegmentInput(Input_Segment, FaultCount, RorT) 

        ############# Flip if needed and get Unit Vectors (triangle Only) ########
        if RorT == "T"
            P1, P2, P3, UnitVector_Normal, UnitVector_StrikeSlip, UnitVector_DipSlip, UnitVector_Slip = 
                RotVerts_UnitVectors(Input_Segment, FaultCount, FaultRakeAngle) 
        end


        if RorT == "R"
            ShearStiffness_H, NormalStiffness_H, Ranks_Shear, Ranks_Normal = 
                HmatBuild_R(ShearModulus, PoissonRatio, ElementRange_SR,Input_Segment, Admissible, Tolerance)  
        elseif RorT == "T"
            ShearStiffness_H, NormalStiffness_H, Ranks_Shear, Ranks_Normal = 
                HmatBuild_T(ShearModulus, PoissonRatio, ElementRange_SR, FaultRakeAngle, FaultCenter, UnitVector_Normal, 
                            UnitVector_Slip, Admissible, Tolerance, P1, P2, P3)   
        end


        save(OutputFileName, 
        "BlockI", BlockI,
        "BlockF", BlockF,
        "Ranks_Shear", Ranks_Shear,
        "Ranks_Normal", Ranks_Normal,
        "ShearStiffness_H",ShearStiffness_H, 
        "NormalStiffness_H", NormalStiffness_H)
        
        println("Saved File Name: ",OutputFileName)


    end
    Discritize()

    """

    fname = "scripts/temp_Discretization/Part$(DistributeIndex).jl"
    open(fname, "w") do file
        write(file, CodeScript )
    end

end




ScriptCombine = """


    using DelimitedFiles
    using Base
    using PyPlot
    using PyCall
    using JLD2
    using LowRankApprox
    using Clustering
    using LinearAlgebra
    @pyimport matplotlib.patches as patches

    pygui(true)

    include("../Functions_BuildInputFile.jl")
    include("../Functions_OKADA3D.jl")
    include("../Functions_Hmatrix.jl")

    OutputFileName="Input_Discretized.jld2"
    InputBulkFileName="Input_BulkFaultGeometry.txt"
    HMatrixStructureFile = "Input_HmatrixStructure.jld2"

    Hmatrix = true
    ShearStiffness_H =[]
    NormalStiffness_H = []
    Ranks_Shear = []
    Ranks_Normal = []

    for Partindex=1:$(HowManyDistribution)

    LoadFileName = "scripts/temp_Discretization/HMatPart_\$(Partindex).jld2"

    BlockI, BlockF, ShearStiffness_H_Part, NormalStiffness_H_Part, Ranks_Shear_Part, Ranks_Normal_Part = 
        load(LoadFileName, "BlockI", "BlockF", "ShearStiffness_H", "NormalStiffness_H",
        "Ranks_Shear", "Ranks_Normal")
        global ShearStiffness_H = [ShearStiffness_H ; ShearStiffness_H_Part]
        global NormalStiffness_H = [NormalStiffness_H ; NormalStiffness_H_Part]
        global Ranks_Shear = [Ranks_Shear ; Ranks_Shear_Part]
        global Ranks_Normal = [Ranks_Normal ; Ranks_Normal_Part]
    end


    ######################## Read General Input File #######################
    LoadingFaultCount, Admissible, ElementRange_SR, Switch_StrikeSlip_or_ReverseNormal, Input_Segment,
    ShearModulus, PoissonRatio, RockDensity,   MinimumNS, RorT, Input_Segment =
        load(HMatrixStructureFile, "LoadingFaultCount", "Admissible", "ElementRange_SR", "Switch_StrikeSlip_or_ReverseNormal", "Input_Segment",
            "ShearModulus", "PoissonRatio", "RockDensity",   "MinimumNS", "RorT", "Input_Segment")
    FaultCount = length(Input_Segment[:,1])

    ######################## Read Segmented Input File #######################
    FaultCenter, Fault_a, Fault_b, Fault_Dc, Fault_Theta_i, Fault_V_i, Fault_Friction_i, Fault_NormalStress, 
    Fault_V_Const, Fault_BulkIndex, FaultLengthStrike, FaultLengthDip, FaultStrikeAngle, 
    FaultDipAngle, FaultRakeAngle, FaultLengthStrike_Bulk, FaultLengthDip_Bulk, NormalStiffnessZero = 
        ReadSegmentInput(Input_Segment, FaultCount, RorT) 

    ############# Flip if needed and get Unit Vectors (triangle Only) ########
    if RorT == "T"
        P1, P2, P3, UnitVector_Normal, UnitVector_StrikeSlip, UnitVector_DipSlip, UnitVector_Slip = 
            RotVerts_UnitVectors(Input_Segment, FaultCount, FaultRakeAngle) 
    end



    save(OutputFileName, 
    "FaultCenter", FaultCenter,
    "ShearModulus", ShearModulus, "RockDensity", RockDensity, "PoissonRatio", PoissonRatio,
    "FaultLengthStrike", FaultLengthStrike, "FaultLengthDip", FaultLengthDip, "FaultStrikeAngle", FaultStrikeAngle, 
    "FaultDipAngle", FaultDipAngle, "FaultRakeAngle", FaultRakeAngle, "Fault_a", Fault_a, "Fault_b", Fault_b, "Fault_Dc", Fault_Dc, 
    "Fault_Theta_i", Fault_Theta_i, "Fault_V_i", Fault_V_i, "Fault_Friction_i", Fault_Friction_i, "Fault_NormalStress", Fault_NormalStress, 
    "Fault_V_Const", Fault_V_Const, "Fault_BulkIndex", Fault_BulkIndex, "FaultLengthStrike_Bulk", FaultLengthStrike_Bulk, 
    "FaultLengthDip_Bulk", FaultLengthDip_Bulk, "FaultCount", FaultCount, "LoadingFaultCount", LoadingFaultCount, 
    "Switch_StrikeSlip_or_ReverseNormal", Switch_StrikeSlip_or_ReverseNormal, "MinimumNormalStress", MinimumNS,
    "NormalStiffnessZero", NormalStiffnessZero, "RorT", RorT,   "Hmatrix", Hmatrix,
    "ShearStiffness_H", ShearStiffness_H, "NormalStiffness_H", NormalStiffness_H, "Admissible", Admissible, 
    "Ranks_Shear", Ranks_Shear, "Ranks_Normal", Ranks_Normal, "ElementRange_SR", ElementRange_SR    
    )
    println("Saved File Name: ",OutputFileName)

    ######################## Save Vertices for Triangles ##############################
    if RorT == "T"        
        file = jldopen(OutputFileName, "a+")
        write(file, "P1", P1) 
        write(file, "P2", P2) 
        write(file, "P3", P3) 
        close(file)
    end



"""

fname = "scripts/temp_Discretization/Combine.jl"
open(fname, "w") do file
    write(file, ScriptCombine )
end





# Block_Ctr_Diam, Block_Range_Level, Input_Segment, LoadingFaultExist, HowManyDivisionEachLevel,
# MinimumElementsToCut, ArrangePoint, Admissible, ElementRange_SR, Switch_StrikeSlip_or_ReverseNormal, 
# ShearModulus, PoissonRatio, RockDensity, DropCrit, DropCritNormalStressMultiplier, MinimumNS, Input_Bulk =
#     load(HMatrixStructureFile, "Block_Ctr_Diam", "Block_Range_Level", "Input_Segment", "LoadingFaultExist", "HowManyDivisionEachLevel",
#             "MinimumElementsToCut", "ArrangePoint", "Admissible", "ElementRange_SR", "Switch_StrikeSlip_or_ReverseNormal", 
#             "ShearModulus", "PoissonRatio", "RockDensity", "DropCrit", "DropCritNormalStressMultiplier", "MinimumNS", "Input_Bulk")


# NormalStiffnessZero = 0    

# ################################################################################
# ################################### Save Files #################################

# FaultCenter = Input_Segment[:,1:3]
# FaultLengthStrike = Input_Segment[:,4]
# FaultLengthDip = Input_Segment[:,5]
# FaultStrikeAngle = Input_Segment[:,6]
# FaultDipAngle = Input_Segment[:,7]
# FaultRakeAngle = Input_Segment[:,8]
# Fault_a = Input_Segment[:,9]
# Fault_b = Input_Segment[:,10]
# Fault_Dc = Input_Segment[:,11]
# Fault_Theta_i = Input_Segment[:,12]
# Fault_V_i = Input_Segment[:,13]
# Fault_Friction_i = Input_Segment[:,14]
# Fault_NormalStress = Input_Segment[:,15]
# Fault_V_Const = Input_Segment[:,16]
# Fault_BulkIndex = Input_Segment[:,17]
# FaultLengthStrike_Bulk = Input_Segment[:,18]
# FaultLengthDip_Bulk = Input_Segment[:,19]
# FaultCount = length(FaultCenter[:,1]) 
# LoadingFaultCount=length(Fault_V_Const[Fault_V_Const.>0])
# SaveOriginalMatrix = 0


# save(OutputFileName, 
# "FaultCenter", FaultCenter,
# "ShearModulus", ShearModulus, "RockDensity", RockDensity, "PoissonRatio", PoissonRatio,
# "FaultLengthStrike", FaultLengthStrike, "FaultLengthDip", FaultLengthDip, "FaultStrikeAngle", FaultStrikeAngle, 
# "FaultDipAngle", FaultDipAngle, "FaultRakeAngle", FaultRakeAngle, "Fault_a", Fault_a, "Fault_b", Fault_b, "Fault_Dc", Fault_Dc, 
# "Fault_Theta_i", Fault_Theta_i, "Fault_V_i", Fault_V_i, "Fault_Friction_i", Fault_Friction_i, "Fault_NormalStress", Fault_NormalStress, 
# "Fault_V_Const", Fault_V_Const, "Fault_BulkIndex", Fault_BulkIndex, "FaultLengthStrike_Bulk", FaultLengthStrike_Bulk, 
# "FaultLengthDip_Bulk", FaultLengthDip_Bulk, "FaultCount", FaultCount, "LoadingFaultCount", LoadingFaultCount, 
# "Switch_StrikeSlip_or_ReverseNormal", Switch_StrikeSlip_or_ReverseNormal, "MinimumNormalStress", MinimumNS,
# "Ranks_Shear", Ranks_Shear, "Ranks_Normal",Ranks_Normal,"ElementRange_SR", ElementRange_SR, "ShearStiffness_H",ShearStiffness_H, "NormalStiffness_H", NormalStiffness_H, "Admissible", Admissible,
# "NormalStiffnessZero", NormalStiffnessZero,"SaveOriginalMatrix",SaveOriginalMatrix)
# println("Saved File Name: ",OutputFileName)


# # println("Saved File Name: ",OutputFileName)
















#=
for DistributeIndex = 1 : HowManyDistribution

    CodeScript = """


    using DelimitedFiles
    using Base
    using PyPlot
    using PyCall
    using JLD2
    using LowRankApprox
    using Clustering
    using LinearAlgebra
    # using Distributed
    @pyimport matplotlib.patches as patches
    
    pygui(true)
    
    include("../Functions_BuildInputFile.jl")
    include("../Functions_OKADA3D.jl")
    include("../Results/Functions_Plot.jl")
    include("../Functions_Hmatrix.jl")
        

    function BuildInputFromBulkGeometry_H()
        BlockI = $(DivisionBlockItoF[DistributeIndex,1])
        BlockF = $(DivisionBlockItoF[DistributeIndex,2])
                    
        InputBulkFileName="../Input_BulkFaultGeometry.txt"
        OutputFileName="HMatPart_$(DistributeIndex).jld2"
        HMatrixStructureFile = "../Input_HmatrixStructure.jld2"



        Block_Ctr_Diam, Block_Range_Level, Input_Segment, LoadingFaultExist, HowManyDivisionEachLevel,
        MinimumElementsToCut, ArrangePoint, Admissible, ElementRange_SR, Switch_StrikeSlip_or_ReverseNormal, 
        ShearModulus, PoissonRatio, RockDensity, DropCrit, DropCritNormalStressMultiplier, MinimumNS, Input_Bulk =
            load(HMatrixStructureFile, "Block_Ctr_Diam", "Block_Range_Level", "Input_Segment", "LoadingFaultExist", "HowManyDivisionEachLevel",
                    "MinimumElementsToCut", "ArrangePoint", "Admissible", "ElementRange_SR", "Switch_StrikeSlip_or_ReverseNormal", 
                    "ShearModulus", "PoissonRatio", "RockDensity", "DropCrit", "DropCritNormalStressMultiplier", "MinimumNS", "Input_Bulk")
    


        
        #####----- Hmatrix compression Tolerance ----#####
        Tolerance = 1e3 # pascal for 1m slip (More approximaion for higher Tolerance). 

        ##########################################################################
        ##########################################################################

        #####---------   Plots?  --------#####
        Plot_HMatStructure = 0 # HMatrix structure plot
        Plot3D_HMatrixBlock = 0 # Blocks in 3D plot
        Plot3D_DiscretizedBlock = 0 # 3D discritized block
        
        ######################################

        ElementPartRoughCount = 2000
        
        if isfile(HMatrixStructureFile)            
            println("H Matrix will be built")
            Block_Ctr_Diam, Block_Range_Level, Input_Segment, LoadingFaultExist, HowManyDivisionEachLevel,
            MinimumElementsToCut, ArrangePoint, Admissible, ElementRange_SR, Switch_StrikeSlip_or_ReverseNormal, 
            ShearModulus, PoissonRatio, RockDensity, DropCrit, DropCritNormalStressMultiplier, MinimumNS, Input_Bulk =
                load(HMatrixStructureFile, "Block_Ctr_Diam", "Block_Range_Level", "Input_Segment", "LoadingFaultExist", "HowManyDivisionEachLevel",
                     "MinimumElementsToCut", "ArrangePoint", "Admissible", "ElementRange_SR", "Switch_StrikeSlip_or_ReverseNormal", 
                     "ShearModulus", "PoissonRatio", "RockDensity", "DropCrit", "DropCritNormalStressMultiplier", "MinimumNS", "Input_Bulk")

            Input_BulkfromFile=readdlm(InputBulkFileName)
            if Input_Bulk != Input_BulkfromFile
                println("???????????????????????????????????????????????????????????????????????????????????????????")
                println("???          H Matrix Structure File dosen't correspond to Bulk Input File              ???")
                println("???????????????????????????????????????????????????????????????????????????????????????????")
            end
            FaultCount = length(Input_Segment[:,1])


        else 
            # println("H Matrix Structure File does not exist")
            println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            println("!!!                        H Matrix Structure File does not exist                       !!!")
            println("!!!      To use Hmatrix apploximation, Run 'RUN_BUILD_HMatrix_Structure.jl' first       !!!")
            println("!!!                          Otherwise, set  HMatrixCompress = 0                        !!!")
            println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            error()

        end



        ####################################################################################
        ############################# Build Stiffness Matrix ###############################
             BlockCount = length(ElementRange_SR[:,1])
            ShearStiffness_H = Any[0]
            NormalStiffness_H = Any[0]
            Ranks_Shear = zeros(Int, BlockCount)
            Ranks_Normal = zeros(Int, BlockCount)

            TotalElments = FaultCount * FaultCount 
            println("Building Hmatrix Block by Block")
            println("Full Matrix Will not be saved")
            println("Preparing for discretization")
            BlockSize = 0.0
            
            # @sync Threads.@threads for BlockIndex = 1: BlockCount
            for BlockIndex = BlockI: BlockF
                
                Input_SegmentS = Input_Segment[ElementRange_SR[BlockIndex,1]:ElementRange_SR[BlockIndex,2],:]
                Input_SegmentR = Input_Segment[ElementRange_SR[BlockIndex,3]:ElementRange_SR[BlockIndex,4],:]
                SourceCount = length(Input_SegmentS[:,1])
                ReceiverCount = length(Input_SegmentR[:,1])
                StiffnessMatrixShearThisBlock=zeros(ReceiverCount,SourceCount)
                StiffnessMatrixNormalThisBlock=zeros(ReceiverCount,SourceCount)

                DivisionCountS = round(Int,SourceCount / ElementPartRoughCount)
                DivisionCountR = round(Int,ReceiverCount / ElementPartRoughCount)
                if DivisionCountS == 0; DivisionCountS =1; end
                if DivisionCountR == 0; DivisionCountR =1; end
                PartedElementCountS = SourceCount ÷ DivisionCountS
                PartedElementCountR = ReceiverCount ÷ DivisionCountR
                TotalParts = DivisionCountS * DivisionCountR
                CurrentPart = 0
                for i=1:DivisionCountS
                    for j=1:DivisionCountR
                        CurrentPart =  CurrentPart +1
                        Init_S = (i-1)*PartedElementCountS + 1
                        Fin_S = i*PartedElementCountS
                        Init_R =  (j-1)*PartedElementCountR + 1
                        Fin_R = j*PartedElementCountR
                        if i == DivisionCountS; Fin_S = SourceCount; end
                        if j == DivisionCountR; Fin_R = ReceiverCount; end
                        BlockSize = BlockSize + (Fin_S - Init_S) * (Fin_R - Init_R)
                        println("Part: ", CurrentPart,"/",DivisionCountS*DivisionCountR, " BlockIndex: ",BlockIndex, "/",BlockCount," Progress:",BlockSize/TotalElments)

                        StiffnessMatrixShearThisBlock[Init_R:Fin_R,Init_S:Fin_S], StiffnessMatrixNormalThisBlock[Init_R:Fin_R,Init_S:Fin_S] = 
                        StiffnessMatrix_ByParts_Calculation_StrikeSlip(Input_SegmentS[Init_S:Fin_S,:], Input_SegmentR[Init_R:Fin_R,:], ShearModulus, PoissonRatio,
                                                            CurrentPart, TotalParts)

                    end
                end                        
                
                ######################################################################
                ############################# Compress ###############################
                if Admissible[BlockIndex] > 0
                    ApproxMatrixS = pqrfact(StiffnessMatrixShearThisBlock, atol = Tolerance)
                    push!(ShearStiffness_H,ApproxMatrixS)
                    Ranks_Shear[BlockIndex] = size(ApproxMatrixS[:Q],2)
                    
                    ApproxMatrixN = pqrfact(StiffnessMatrixNormalThisBlock, atol = Tolerance)
                    push!(NormalStiffness_H,ApproxMatrixN)
                    Ranks_Normal[BlockIndex] = size(ApproxMatrixN[:Q],2)
                else 
                    push!(ShearStiffness_H,StiffnessMatrixShearThisBlock)
                    push!(NormalStiffness_H,StiffnessMatrixNormalThisBlock)
                end
            end
            ShearStiffness_H = ShearStiffness_H[2:end]
            NormalStiffness_H = NormalStiffness_H[2:end]
            


        save(OutputFileName, 
        "BlockI", BlockI,
        "BlockF", BlockF,
        "Ranks_Shear", Ranks_Shear,
        "Ranks_Normal", Ranks_Normal,
        "ShearStiffness_H",ShearStiffness_H, 
        "NormalStiffness_H", NormalStiffness_H)
        
        println("Saved File Name: ",OutputFileName)



    end


    Input = BuildInputFromBulkGeometry_H()

    """

    fname = "Part$(DistributeIndex).jl"
    open(fname, "w") do file
        write(file, CodeScript )
    end


end



ScriptCombine = """


using DelimitedFiles
using Base
using PyPlot
using PyCall
using JLD2
using LowRankApprox
using Clustering
using LinearAlgebra
using Distributed
@pyimport matplotlib.patches as patches

pygui(true)

include("../Functions_BuildInputFile.jl")
include("../Functions_OKADA3D.jl")
include("../Results/Functions_Plot.jl")
include("../Functions_Hmatrix.jl")

OutputFileName="Input_Discretized.jld2"
InputBulkFileName="../Input_BulkFaultGeometry.txt"
HMatrixStructureFile = "../Input_HmatrixStructure.jld2"





ShearStiffness_H =[]
NormalStiffness_H = []
Ranks_Shear = zeros($(TotalBlock))
Ranks_Normal = zeros($(TotalBlock))

for Partindex=1:$(HowManyDistribution)

LoadFileName = "HMatPart_\$(Partindex).jld2"

BlockI, BlockF, ShearStiffness_H_Part, NormalStiffness_H_Part, Ranks_Shear_Part, Ranks_Normal_Part = 
    load(LoadFileName, "BlockI", "BlockF", "ShearStiffness_H", "NormalStiffness_H",
    "Ranks_Shear", "Ranks_Normal")
    global ShearStiffness_H = [ShearStiffness_H ; ShearStiffness_H_Part]
    global NormalStiffness_H = [NormalStiffness_H ; NormalStiffness_H_Part]
    global Ranks_Shear = Ranks_Shear .+ Ranks_Shear_Part
    global Ranks_Normal = Ranks_Normal .+ Ranks_Normal_Part
end



Block_Ctr_Diam, Block_Range_Level, Input_Segment, LoadingFaultExist, HowManyDivisionEachLevel,
MinimumElementsToCut, ArrangePoint, Admissible, ElementRange_SR, Switch_StrikeSlip_or_ReverseNormal, 
ShearModulus, PoissonRatio, RockDensity, DropCrit, DropCritNormalStressMultiplier, MinimumNS, Input_Bulk =
    load(HMatrixStructureFile, "Block_Ctr_Diam", "Block_Range_Level", "Input_Segment", "LoadingFaultExist", "HowManyDivisionEachLevel",
            "MinimumElementsToCut", "ArrangePoint", "Admissible", "ElementRange_SR", "Switch_StrikeSlip_or_ReverseNormal", 
            "ShearModulus", "PoissonRatio", "RockDensity", "DropCrit", "DropCritNormalStressMultiplier", "MinimumNS", "Input_Bulk")


NormalStiffnessZero = 0    

################################################################################
################################### Save Files #################################

FaultCenter = Input_Segment[:,1:3]
FaultLengthStrike = Input_Segment[:,4]
FaultLengthDip = Input_Segment[:,5]
FaultStrikeAngle = Input_Segment[:,6]
FaultDipAngle = Input_Segment[:,7]
FaultRakeAngle = Input_Segment[:,8]
Fault_a = Input_Segment[:,9]
Fault_b = Input_Segment[:,10]
Fault_Dc = Input_Segment[:,11]
Fault_Theta_i = Input_Segment[:,12]
Fault_V_i = Input_Segment[:,13]
Fault_Friction_i = Input_Segment[:,14]
Fault_NormalStress = Input_Segment[:,15]
Fault_V_Const = Input_Segment[:,16]
Fault_BulkIndex = Input_Segment[:,17]
FaultLengthStrike_Bulk = Input_Segment[:,18]
FaultLengthDip_Bulk = Input_Segment[:,19]
FaultCount = length(FaultCenter[:,1]) 
LoadingFaultCount=length(Fault_V_Const[Fault_V_Const.>0])
SaveOriginalMatrix = 0


save(OutputFileName, 
"FaultCenter", FaultCenter,
"ShearModulus", ShearModulus, "RockDensity", RockDensity, "PoissonRatio", PoissonRatio,
"FaultLengthStrike", FaultLengthStrike, "FaultLengthDip", FaultLengthDip, "FaultStrikeAngle", FaultStrikeAngle, 
"FaultDipAngle", FaultDipAngle, "FaultRakeAngle", FaultRakeAngle, "Fault_a", Fault_a, "Fault_b", Fault_b, "Fault_Dc", Fault_Dc, 
"Fault_Theta_i", Fault_Theta_i, "Fault_V_i", Fault_V_i, "Fault_Friction_i", Fault_Friction_i, "Fault_NormalStress", Fault_NormalStress, 
"Fault_V_Const", Fault_V_Const, "Fault_BulkIndex", Fault_BulkIndex, "FaultLengthStrike_Bulk", FaultLengthStrike_Bulk, 
"FaultLengthDip_Bulk", FaultLengthDip_Bulk, "FaultCount", FaultCount, "LoadingFaultCount", LoadingFaultCount, 
"Switch_StrikeSlip_or_ReverseNormal", Switch_StrikeSlip_or_ReverseNormal, "MinimumNormalStress", MinimumNS,
"Ranks_Shear", Ranks_Shear, "Ranks_Normal",Ranks_Normal,"ElementRange_SR", ElementRange_SR, "ShearStiffness_H",ShearStiffness_H, "NormalStiffness_H", NormalStiffness_H, "Admissible", Admissible,
"NormalStiffnessZero", NormalStiffnessZero,"SaveOriginalMatrix",SaveOriginalMatrix)
println("Saved File Name: ",OutputFileName)


# println("Saved File Name: ",OutputFileName)

"""


fname = "Combine.jl"
open(fname, "w") do file
    write(file, ScriptCombine )
end



#=


function BuildInputFromBulkGeometry_H()

    InputBulkFileName="../Input_BulkFaultGeometry.txt"
    OutputFileName="Input_Discretized.jld2"
    HMatrixStructureFile = "../Input_HmatrixStructure.jld2"
    ##########################################################################
    ########################## Hmatrix compress? #############################
    HMatrixCompress = 1 # If this is 1, stiffness Matrix will be compressed using Input_HmatrixStructure.jld2
    SaveOriginalMatrix = 1  # 1: save Original Matrix (can be very large), 0: Discard Original Matrix. 
    
    #####----- Hmatrix compression Tolerance ----#####
    Tolerance = 1e3 # pascal for 1m slip (More approximaion for higher Tolerance). 

    ##########################################################################
    ##########################################################################

    #####---------   Plots?  --------#####
    Plot_HMatStructure = 0 # HMatrix structure plot
    Plot3D_HMatrixBlock = 0 # Blocks in 3D plot
    Plot3D_DiscretizedBlock = 0 # 3D discritized block
    
    ######################################

    ElementPartRoughCount = 2000
    
    if isfile(HMatrixStructureFile)            
        println("H Matrix will be built")
        Block_Ctr_Diam, Block_Range_Level, Input_Segment, LoadingFaultExist, HowManyDivisionEachLevel,
        MinimumElementsToCut, ArrangePoint, Admissible, ElementRange_SR, Switch_StrikeSlip_or_ReverseNormal, 
        ShearModulus, PoissonRatio, RockDensity, DropCrit, DropCritNormalStressMultiplier, MinimumNS, Input_Bulk =
            load(HMatrixStructureFile, "Block_Ctr_Diam", "Block_Range_Level", "Input_Segment", "LoadingFaultExist", "HowManyDivisionEachLevel",
                    "MinimumElementsToCut", "ArrangePoint", "Admissible", "ElementRange_SR", "Switch_StrikeSlip_or_ReverseNormal", 
                    "ShearModulus", "PoissonRatio", "RockDensity", "DropCrit", "DropCritNormalStressMultiplier", "MinimumNS", "Input_Bulk")

        Input_BulkfromFile=readdlm(InputBulkFileName)
        if Input_Bulk != Input_BulkfromFile
            println("???????????????????????????????????????????????????????????????????????????????????????????")
            println("???          H Matrix Structure File dosen't correspond to Bulk Input File              ???")
            println("???????????????????????????????????????????????????????????????????????????????????????????")
        end
        FaultCount = length(Input_Segment[:,1])


    else 
        # println("H Matrix Structure File does not exist")
        println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        println("!!!                        H Matrix Structure File does not exist                       !!!")
        println("!!!      To use Hmatrix apploximation, Run 'RUN_BUILD_HMatrix_Structure.jl' first       !!!")
        println("!!!                          Otherwise, set  HMatrixCompress = 0                        !!!")
        println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        error()

    end



    ####################################################################################
    ############################# Build Stiffness Matrix ###############################
    BlockCount = length(ElementRange_SR[:,1])
    ShearStiffness_H = Any[0]
    NormalStiffness_H = Any[0]
    Ranks_Shear = zeros(Int, BlockCount)
    Ranks_Normal = zeros(Int, BlockCount)

    TotalElments = FaultCount * FaultCount 
    println("Building Hmatrix Block by Block")
    println("Full Matrix Will not be saved")
    println("Preparing for discretization")
    BlockSize = 0.0
    
    # @sync Threads.@threads for BlockIndex = 1: BlockCount
    @time for BlockIndex = 1: BlockCount
        
        Input_SegmentS = Input_Segment[ElementRange_SR[BlockIndex,1]:ElementRange_SR[BlockIndex,2],:]
        Input_SegmentR = Input_Segment[ElementRange_SR[BlockIndex,3]:ElementRange_SR[BlockIndex,4],:]
        SourceCount = length(Input_SegmentS[:,1])
        ReceiverCount = length(Input_SegmentR[:,1])
        StiffnessMatrixShearThisBlock=zeros(ReceiverCount,SourceCount)
        StiffnessMatrixNormalThisBlock=zeros(ReceiverCount,SourceCount)

        DivisionCountS = round(Int,SourceCount / ElementPartRoughCount)
        DivisionCountR = round(Int,ReceiverCount / ElementPartRoughCount)
        if DivisionCountS == 0; DivisionCountS =1; end
        if DivisionCountR == 0; DivisionCountR =1; end
        PartedElementCountS = SourceCount ÷ DivisionCountS
        PartedElementCountR = ReceiverCount ÷ DivisionCountR
        TotalParts = DivisionCountS * DivisionCountR
        CurrentPart = 0
        for i=1:DivisionCountS
            for j=1:DivisionCountR
                CurrentPart =  CurrentPart +1
                Init_S = (i-1)*PartedElementCountS + 1
                Fin_S = i*PartedElementCountS
                Init_R =  (j-1)*PartedElementCountR + 1
                Fin_R = j*PartedElementCountR
                if i == DivisionCountS; Fin_S = SourceCount; end
                if j == DivisionCountR; Fin_R = ReceiverCount; end
                BlockSize = BlockSize + (Fin_S - Init_S) * (Fin_R - Init_R)
                println("Part: ", CurrentPart,"/",DivisionCountS*DivisionCountR, " BlockIndex: ",BlockIndex, "/",BlockCount," Progress:",BlockSize/TotalElments)

                if Switch_StrikeSlip_or_ReverseNormal == 1
                StiffnessMatrixShearThisBlock[Init_R:Fin_R,Init_S:Fin_S], StiffnessMatrixNormalThisBlock[Init_R:Fin_R,Init_S:Fin_S] = 
                StiffnessMatrix_ByParts_Calculation_StrikeSlip(Input_SegmentS[Init_S:Fin_S,:], Input_SegmentR[Init_R:Fin_R,:], ShearModulus, PoissonRatio,
                                                    CurrentPart, TotalParts)

                elseif Switch_StrikeSlip_or_ReverseNormal == 2                                                        
                StiffnessMatrixShearThisBlock[Init_R:Fin_R,Init_S:Fin_S], StiffnessMatrixNormalThisBlock[Init_R:Fin_R,Init_S:Fin_S] = 
                StiffnessMatrix_ByParts_Calculation_NormalReverse(Input_SegmentS[Init_S:Fin_S,:], Input_SegmentR[Init_R:Fin_R,:], ShearModulus, PoissonRatio,
                                                    CurrentPart, TotalParts)                                                                    
                else
                    println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                    println("!!!!!!! Switch_StrikeSlip_or_ReverseNormal should be either 1 or 2  !!!!!!!")
                    println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                end
            end
        end                        
        
        ######################################################################
        ############################# Compress ###############################
        if Admissible[BlockIndex] > 0
            ApproxMatrixS = pqrfact(StiffnessMatrixShearThisBlock, atol = Tolerance)
            push!(ShearStiffness_H,ApproxMatrixS)
            Ranks_Shear[BlockIndex] = size(ApproxMatrixS[:Q],2)
            
            ApproxMatrixN = pqrfact(StiffnessMatrixNormalThisBlock, atol = Tolerance)
            push!(NormalStiffness_H,ApproxMatrixN)
            Ranks_Normal[BlockIndex] = size(ApproxMatrixN[:Q],2)
        else 
            push!(ShearStiffness_H,StiffnessMatrixShearThisBlock)
            push!(NormalStiffness_H,StiffnessMatrixNormalThisBlock)
        end
    end
    ShearStiffness_H = ShearStiffness_H[2:end]
    NormalStiffness_H = NormalStiffness_H[2:end]



    NormalStiffnessZero = 0    

    ################################################################################
    ################################### Save Files #################################

    FaultCenter = Input_Segment[:,1:3]
    FaultLengthStrike = Input_Segment[:,4]
    FaultLengthDip = Input_Segment[:,5]
    FaultStrikeAngle = Input_Segment[:,6]
    FaultDipAngle = Input_Segment[:,7]
    FaultLLRR = Input_Segment[:,8]
    Fault_a = Input_Segment[:,9]
    Fault_b = Input_Segment[:,10]
    Fault_Dc = Input_Segment[:,11]
    Fault_Theta_i = Input_Segment[:,12]
    Fault_V_i = Input_Segment[:,13]
    Fault_Friction_i = Input_Segment[:,14]
    Fault_NormalStress = Input_Segment[:,15]
    Fault_V_Const = Input_Segment[:,16]
    Fault_BulkIndex = Input_Segment[:,17]
    FaultLengthStrike_Bulk = Input_Segment[:,18]
    FaultLengthDip_Bulk = Input_Segment[:,19]

    LoadingFaultCount=length(Fault_V_Const[Fault_V_Const.>0])

            save(OutputFileName, 
            "FaultCenter", FaultCenter,
            "BlockI", BlockI,
            "BlockF", BlockF,
            "ShearStiffness_H",ShearStiffness_H, 
            "NormalStiffness_H", NormalStiffness_H)
            println("Saved File Name: ",OutputFileName)




end


Input = BuildInputFromBulkGeometry_H()


=#
=#