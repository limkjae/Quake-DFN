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

include("Functions_BuildInputFile.jl")
include("Functions_OKADA3D.jl")
include("Results/Functions_Plot.jl")
include("Functions_Hmatrix.jl")


AAA = @isdefined FunctionRead
if AAA == false
    include("TriangleDev/Functions_TDstressHS.jl")
    FunctionRead = 1
end

InputBulkFileName="Input_BulkFaultGeometry.txt"
OutputFileName="Input_Discretized.jld2"
HMatrixStructureFile = "Input_HmatrixStructure.jld2"
##########################################################################
########################## Hmatrix compress? #############################
HMatrixCompress = 1 # If this is 1, stiffness Matrix will be compressed using Input_HmatrixStructure.jld2
SaveOriginalMatrix = 0  # 1: save Original Matrix (can be very large), 0: Discard Original Matrix. 

#####----- Hmatrix compression Tolerance ----#####
Tolerance = 1e3 # pascal for 1m slip (More approximaion for higher Tolerance). 

##########################################################################
##########################################################################

#####---------   Plots?  --------#####
Plot_HMatStructure = 1 # HMatrix structure plot
Plot3D_HMatrixBlock = 0 # Blocks in 3D plot
Plot3D_DiscretizedBlock = 0 # 3D discritized block

######################################



######################## Read Bulk Input File ########################
FaultCenter, Fault_a, Fault_b, Fault_Dc, Fault_Theta_i, Fault_V_i, Fault_Friction_i, Fault_NormalStress, 
Fault_V_Const, Fault_BulkIndex, FaultLengthStrike, FaultLengthDip, FaultStrikeAngle, 
FaultDipAngle, FaultRakeAngle, FaultLengthStrike_Bulk, FaultLengthDip_Bulk, NormalStiffnessZero,
Input_Segment, LoadingFaultCount, ShearModulus, PoissonRatio, RockDensity, 
Switch_StrikeSlip_or_ReverseNormal, DropCrit, DropCritNormalStressMultiplier, MinimumNS, RorT, FaultCount =
    ReadBulkInput(InputBulkFileName)




######################## Flip if needed and get Unit Vectors ########################

P1, P2, P3, UnitVector_Normal, UnitVector_StrikeSlip, UnitVector_DipSlip, UnitVector_Slip = 
    RotVerts_UnitVectors(Input_Segment, FaultCount, FaultRakeAngle) 


###################### Build Stiffness Matrix ########################

ElementPartRoughCount = 2000
println("building full stiffness matrix (No-Hmatrix Compression)")

StiffnessMatrix_Shear, StiffnessMatrix_Normal = 
    BuildMatrixByParts_Even_Tri(P1, P2, P3, FaultRakeAngle, FaultCenter, UnitVector_Normal, UnitVector_Slip,
        FaultCount, ElementPartRoughCount, ShearModulus, PoissonRatio) 

println("Stiffness Matrix Built")


save(OutputFileName, 
"StiffnessMatrixShear", StiffnessMatrix_Shear, "StiffnessMatrixNormal", StiffnessMatrix_Normal, 
"FaultCenter", FaultCenter, "P1", P1, "P2", P2, "P3", P3,
"ShearModulus", ShearModulus, "RockDensity", RockDensity, "PoissonRatio", PoissonRatio,
"FaultLengthStrike", FaultLengthStrike, "FaultLengthDip", FaultLengthDip, "FaultStrikeAngle", FaultStrikeAngle, 
"FaultDipAngle", FaultDipAngle, "FaultRakeAngle", FaultRakeAngle, "Fault_a", Fault_a, "Fault_b", Fault_b, "Fault_Dc", Fault_Dc, 
"Fault_Theta_i", Fault_Theta_i, "Fault_V_i", Fault_V_i, "Fault_Friction_i", Fault_Friction_i, "Fault_NormalStress", Fault_NormalStress, 
"Fault_V_Const", Fault_V_Const, "Fault_BulkIndex", Fault_BulkIndex, "FaultLengthStrike_Bulk", FaultLengthStrike_Bulk, 
"FaultLengthDip_Bulk", FaultLengthDip_Bulk, "FaultCount", FaultCount, "LoadingFaultCount", LoadingFaultCount,
"Switch_StrikeSlip_or_ReverseNormal", Switch_StrikeSlip_or_ReverseNormal, "MinimumNormalStress", MinimumNS,
"NormalStiffnessZero", NormalStiffnessZero)
println("Saved File Name: ",OutputFileName)

