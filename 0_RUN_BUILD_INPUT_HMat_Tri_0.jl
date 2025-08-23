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

AAA = @isdefined FunctionRead
if AAA == false
    include("TriangleDev/Functions_TDstressHS.jl")
    FunctionRead = 1
end

include("Functions_BuildInputFile.jl")
include("Functions_OKADA3D.jl")
include("Results/Functions_Plot.jl")
include("Functions_Hmatrix.jl")

InputBulkFileName="Input_BulkFaultGeometry.txt"
OutputFileName="Input_Discretized.jld2"
HMatrixStructureFile = "Input_HmatrixStructure.jld2"



##########################################################################
##########################       Input       #############################

Hmatrix = true # true: build Hmatrix structure from Input_HmatrixStructure.jld2, false: build full matrix without compression

Tolerance = 1e3 #  Hmatrix compression Tolerance. pascal for 1m slip (More approximaion for higher Tolerance). 

#####---------   Plots?  --------#####
Plot_HMatStructure = 1 # HMatrix structure plot
Plot3D_HMatrixBlock = 0 # Blocks in 3D plot
Plot3D_DiscretizedBlock = 0 # 3D discritized block

##########################################################################
##########################################################################



######################## Read Bulk Input File ############################
if Hmatrix == false

    Input_Segment, LoadingFaultCount, ShearModulus, PoissonRatio, RockDensity, 
    Switch_StrikeSlip_or_ReverseNormal, DropCrit, DropCritNormalStressMultiplier, MinimumNS, RorT, FaultCount =
        ReadBulkInput(InputBulkFileName)

elseif Hmatrix == true
    println("H Matrix will be built")
    Block_Ctr_Diam, Block_Range_Level, Input_Segment, LoadingFaultExist, LoadingFaultCount,
    MinimumElementsToCut, ArrangePoint, Admissible, ElementRange_SR, Switch_StrikeSlip_or_ReverseNormal, 
    ShearModulus, PoissonRatio, RockDensity, DropCrit, DropCritNormalStressMultiplier, MinimumNS, RorT =
        load(HMatrixStructureFile, "Block_Ctr_Diam", "Block_Range_Level", "Input_Segment", "LoadingFaultExist", "LoadingFaultCount",
            "MinimumElementsToCut", "ArrangePoint", "Admissible", "ElementRange_SR", "Switch_StrikeSlip_or_ReverseNormal", 
            "ShearModulus", "PoissonRatio", "RockDensity", "DropCrit", "DropCritNormalStressMultiplier", "MinimumNS", "RorT")

    FaultCount = length(Input_Segment[:,1])

end



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



################################ HMat Discritize #############################
            println("Build Full Matrix and then compress")

            StiffnessMatrixShearOriginal=zeros(FaultCount,FaultCount)
            StiffnessMatrixNormalOriginal=zeros(FaultCount,FaultCount)    

            println("Preparing for discretization")
            StiffnessMatrixShearOriginal, StiffnessMatrixNormalOriginal = 
                BuildMatrixByParts_Even_Tri(P1, P2, P3, FaultRakeAngle, FaultCenter, UnitVector_Normal, UnitVector_Slip,
                    FaultCount, ElementPartRoughCount, ShearModulus, PoissonRatio) 
                    
            BlockCount = length(ElementRange_SR[:,1])
            ShearStiffness_H = Any[0]
            Ranks_Shear = zeros(Int, BlockCount)
            NormalStiffness_H = Any[0]
            Ranks_Normal = zeros(Int, BlockCount)
            BlockIndex = 0
            println("compressing")
            for i=1:BlockCount
                BlockIndex = BlockIndex + 1
                
                if Admissible[BlockIndex] > 0
                    OrigianlMatrixToApproximate = StiffnessMatrixShearOriginal[ElementRange_SR[i,3]:ElementRange_SR[i,4],ElementRange_SR[i,1]:ElementRange_SR[i,2]]
                    ApproxMatrixS = pqrfact(OrigianlMatrixToApproximate, atol = Tolerance)
                    push!(ShearStiffness_H,ApproxMatrixS)
                    Ranks_Shear[BlockIndex] = size(ApproxMatrixS[:Q],2)
                    
                    OrigianlMatrixToApproximate = StiffnessMatrixNormalOriginal[ElementRange_SR[i,3]:ElementRange_SR[i,4],ElementRange_SR[i,1]:ElementRange_SR[i,2]]
                    ApproxMatrixN = pqrfact(OrigianlMatrixToApproximate, atol = Tolerance)
                    push!(NormalStiffness_H,ApproxMatrixN)
                    Ranks_Normal[BlockIndex] = size(ApproxMatrixN[:Q],2)
                else 
                    push!(ShearStiffness_H,StiffnessMatrixShearOriginal[ElementRange_SR[i,3]:ElementRange_SR[i,4],ElementRange_SR[i,1]:ElementRange_SR[i,2]])
                    push!(NormalStiffness_H,StiffnessMatrixNormalOriginal[ElementRange_SR[i,3]:ElementRange_SR[i,4],ElementRange_SR[i,1]:ElementRange_SR[i,2]])
                end
    
            end
            ShearStiffness_H = ShearStiffness_H[2:end]
            NormalStiffness_H = NormalStiffness_H[2:end]





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
"NormalStiffnessZero", NormalStiffnessZero)
println("Saved File Name: ",OutputFileName)

######################## Save Input Files ##############################
if RorT == "T"        
    file = jldopen(OutputFileName, "a+")
    write(file, "P1", P1) 
    write(file, "P2", P2) 
    write(file, "P3", P3) 
    close(file)

end







# ################################  Discritize ##############################
# ElementPartRoughCount = 2000

# lambda = 2 * ShearModulus * PoissonRatio / (1 - 2 * PoissonRatio)
# BlockCount = length(ElementRange_SR[:,1])
# ShearStiffness_H = Any[0]
# NormalStiffness_H = Any[0]
# Ranks_Shear = zeros(Int, BlockCount)
# Ranks_Normal = zeros(Int, BlockCount)

# TotalElments = FaultCount * FaultCount 
# println("Building Hmatrix Block by Block")
# println("Full Matrix Will not be saved")
# println("Preparing for discretization")
# BlockSize = 0.0

# # @sync Threads.@threads for BlockIndex = 1: BlockCount
# for BlockIndex = 1: BlockCount
#     # for BlockIndex = 1: 5
# # BlockIndex =2
#     Input_SegmentS = Input_Segment[ElementRange_SR[BlockIndex,1]:ElementRange_SR[BlockIndex,2],:]
#     Input_SegmentR = Input_Segment[ElementRange_SR[BlockIndex,3]:ElementRange_SR[BlockIndex,4],:]
#     SourceCount = length(Input_SegmentS[:,1])
#     ReceiverCount = length(Input_SegmentR[:,1])
#     StiffnessMatrixShearThisBlock=zeros(ReceiverCount,SourceCount)
#     StiffnessMatrixNormalThisBlock=zeros(ReceiverCount,SourceCount)

#     DivisionCountS = round(Int,SourceCount / ElementPartRoughCount)
#     DivisionCountR = round(Int,ReceiverCount / ElementPartRoughCount)

#     if DivisionCountS == 0; DivisionCountS =1; end
#     if DivisionCountR == 0; DivisionCountR =1; end
#     PartedElementCountS = SourceCount ÷ DivisionCountS
#     PartedElementCountR = ReceiverCount ÷ DivisionCountR
#     TotalParts = DivisionCountS * DivisionCountR
#     CurrentPart = 0
#     for i=1:DivisionCountS
#         for j=1:DivisionCountR
#             CurrentPart =  CurrentPart +1
#             Init_S = (i-1)*PartedElementCountS + 1
#             Fin_S = i*PartedElementCountS
#             Init_R =  (j-1)*PartedElementCountR + 1
#             Fin_R = j*PartedElementCountR
#             if i == DivisionCountS; Fin_S = SourceCount; end
#             if j == DivisionCountR; Fin_R = ReceiverCount; end
#             BlockSize = BlockSize + (Fin_S - Init_S) * (Fin_R - Init_R)
#             println("Part: ", CurrentPart,"/",DivisionCountS*DivisionCountR, " BlockIndex: ",BlockIndex, "/",BlockCount," Progress:",BlockSize/TotalElments)

#             if RorT == "R"
#                         StiffnessMatrixShearThisBlock[Init_R:Fin_R,Init_S:Fin_S], StiffnessMatrixNormalThisBlock[Init_R:Fin_R,Init_S:Fin_S] = 
#                         StiffnessMatrix_ByParts_Calculation_Rec(Input_SegmentS[Init_S:Fin_S,:], Input_SegmentR[Init_R:Fin_R,:], ShearModulus, PoissonRatio,
#                                                             CurrentPart, TotalParts)
#             else 
#                         StiffnessMatrixShearThisBlock[Init_R:Fin_R,Init_S:Fin_S], StiffnessMatrixNormalThisBlock[Init_R:Fin_R,Init_S:Fin_S] = 
#                         StiffnessMatrix_ByParts_Calculation_Tri(P1[Init_S:Fin_S,:], P2[Init_S:Fin_S,:], P3[Init_S:Fin_S,:], FaultRakeAngle[Init_S:Fin_S,:],
#                                                                 FaultCenter[Init_R:Fin_R,:], ShearModulus, lambda,  
#                                                                 UnitVector_Normal[Init_R:Fin_R,:], UnitVector_Slip[Init_R:Fin_R,:])   
#             end                                               
            
                                            
#         end
#     end                        

#         ######################################################################
#         ############################# Compress ###############################
#     if Admissible[BlockIndex] > 0
#         ApproxMatrixS = pqrfact(StiffnessMatrixShearThisBlock, atol = Tolerance)
#         push!(ShearStiffness_H,ApproxMatrixS)
#         Ranks_Shear[BlockIndex] = size(ApproxMatrixS[:Q],2)
        
#         ApproxMatrixN = pqrfact(StiffnessMatrixNormalThisBlock, atol = Tolerance)
#         push!(NormalStiffness_H,ApproxMatrixN)
#         Ranks_Normal[BlockIndex] = size(ApproxMatrixN[:Q],2)
#     else 
#         push!(ShearStiffness_H,StiffnessMatrixShearThisBlock)
#         push!(NormalStiffness_H,StiffnessMatrixNormalThisBlock)
#     end
# # ApproxMatrixS*ones(200)
# end
#             ShearStiffness_H = ShearStiffness_H[2:end]
#             NormalStiffness_H = NormalStiffness_H[2:end]

# # ShearStiffness_H[2]









# ################################ HMat Discritize #############################
#             println("Build Full Matrix and then compress")

#             StiffnessMatrixShearOriginal=zeros(FaultCount,FaultCount)
#             StiffnessMatrixNormalOriginal=zeros(FaultCount,FaultCount)    

#             println("Preparing for discretization")
#             StiffnessMatrixShearOriginal, StiffnessMatrixNormalOriginal = 
#                 BuildMatrixByParts_Even_Tri(P1, P2, P3, FaultRakeAngle, FaultCenter, UnitVector_Normal, UnitVector_Slip,
#                     FaultCount, ElementPartRoughCount, ShearModulus, PoissonRatio) 
                    
#             BlockCount = length(ElementRange_SR[:,1])
#             ShearStiffness_H = Any[0]
#             Ranks_Shear = zeros(Int, BlockCount)
#             NormalStiffness_H = Any[0]
#             Ranks_Normal = zeros(Int, BlockCount)
#             BlockIndex = 0
#             println("compressing")
#             for i=1:BlockCount
#                 BlockIndex = BlockIndex + 1
                
#                 if Admissible[BlockIndex] > 0
#                     OrigianlMatrixToApproximate = StiffnessMatrixShearOriginal[ElementRange_SR[i,3]:ElementRange_SR[i,4],ElementRange_SR[i,1]:ElementRange_SR[i,2]]
#                     ApproxMatrixS = pqrfact(OrigianlMatrixToApproximate, atol = Tolerance)
#                     push!(ShearStiffness_H,ApproxMatrixS)
#                     Ranks_Shear[BlockIndex] = size(ApproxMatrixS[:Q],2)
                    
#                     OrigianlMatrixToApproximate = StiffnessMatrixNormalOriginal[ElementRange_SR[i,3]:ElementRange_SR[i,4],ElementRange_SR[i,1]:ElementRange_SR[i,2]]
#                     ApproxMatrixN = pqrfact(OrigianlMatrixToApproximate, atol = Tolerance)
#                     push!(NormalStiffness_H,ApproxMatrixN)
#                     Ranks_Normal[BlockIndex] = size(ApproxMatrixN[:Q],2)
#                 else 
#                     push!(ShearStiffness_H,StiffnessMatrixShearOriginal[ElementRange_SR[i,3]:ElementRange_SR[i,4],ElementRange_SR[i,1]:ElementRange_SR[i,2]])
#                     push!(NormalStiffness_H,StiffnessMatrixNormalOriginal[ElementRange_SR[i,3]:ElementRange_SR[i,4],ElementRange_SR[i,1]:ElementRange_SR[i,2]])
#                 end
    
#             end
#             ShearStiffness_H = ShearStiffness_H[2:end]
#             NormalStiffness_H = NormalStiffness_H[2:end]

























# ######################## Build Stiffness Matrix ##########################

# println("building full stiffness matrix (No-Hmatrix Compression)")

# if RorT == "R"
#     println("Preparing for discretization")
#         StiffnessMatrix_Shear, StiffnessMatrix_Normal = 
#             BuildMatrixByParts_Even_Rec(FaultCount, ElementPartRoughCount, Input_Segment,  ShearModulus, PoissonRatio)
            
# else 

#     StiffnessMatrix_Shear, StiffnessMatrix_Normal = 
#         BuildMatrixByParts_Even_Tri(P1, P2, P3, FaultRakeAngle, FaultCenter, UnitVector_Normal, UnitVector_Slip,
#             FaultCount, ElementPartRoughCount, ShearModulus, PoissonRatio) 
# end
# println("Stiffness Matrix Built")

# save(OutputFileName, 
# "StiffnessMatrixShear", StiffnessMatrix_Shear, "StiffnessMatrixNormal", StiffnessMatrix_Normal, "FaultCenter", FaultCenter,
# "ShearModulus", ShearModulus, "RockDensity", RockDensity, "PoissonRatio", PoissonRatio,
# "FaultLengthStrike", FaultLengthStrike, "FaultLengthDip", FaultLengthDip, "FaultStrikeAngle", FaultStrikeAngle, 
# "FaultDipAngle", FaultDipAngle, "FaultRakeAngle", FaultRakeAngle, "Fault_a", Fault_a, "Fault_b", Fault_b, "Fault_Dc", Fault_Dc, 
# "Fault_Theta_i", Fault_Theta_i, "Fault_V_i", Fault_V_i, "Fault_Friction_i", Fault_Friction_i, "Fault_NormalStress", Fault_NormalStress, 
# "Fault_V_Const", Fault_V_Const, "Fault_BulkIndex", Fault_BulkIndex, "FaultLengthStrike_Bulk", FaultLengthStrike_Bulk, 
# "FaultLengthDip_Bulk", FaultLengthDip_Bulk, "FaultCount", FaultCount, "LoadingFaultCount", LoadingFaultCount,
# "Switch_StrikeSlip_or_ReverseNormal", Switch_StrikeSlip_or_ReverseNormal, "MinimumNormalStress", MinimumNS, "RorT", RorT,
# "NormalStiffnessZero", NormalStiffnessZero)
# println("Saved File Name: ",OutputFileName)
