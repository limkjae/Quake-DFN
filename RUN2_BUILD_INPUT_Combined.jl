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
Hmatrix = false # true: build Hmatrix structure from Input_HmatrixStructure.jld2, false: build full matrix without compression
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


if Hmatrix == true
    if RorT == "R"
        ShearStiffness_H, NormalStiffness_H, Ranks_Shear, Ranks_Normal = 
            HmatBuild_R(ShearModulus, PoissonRatio, ElementRange_SR,Input_Segment) 
        
        
    elseif RorT == "T"
        ShearStiffness_H, NormalStiffness_H, Ranks_Shear, Ranks_Normal = 
            HmatBuild_T(ShearModulus, PoissonRatio, ElementRange_SR, FaultRakeAngle, FaultCenter, UnitVector_Normal, UnitVector_Slip) 
    end

elseif Hmatrix == false
    if RorT == "R"
        StiffnessMatrixShear, StiffnessMatrixNormal = 
            BuildMatrixByParts_Even_Rec(FaultCount, Input_Segment,  ShearModulus, PoissonRatio)
        
        
    elseif RorT == "T"
        StiffnessMatrixShear, StiffnessMatrixNormal = 
            BuildMatrixByParts_Even_Tri(P1, P2, P3, FaultRakeAngle, FaultCenter, UnitVector_Normal, UnitVector_Slip,
                                FaultCount, ShearModulus, PoissonRatio)                            
    end
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
"NormalStiffnessZero", NormalStiffnessZero, "RorT", RorT,   "Hmatrix", Hmatrix)
println("Saved File Name: ",OutputFileName)

file = jldopen(OutputFileName, "a+")
if Hmatrix == true
    write(file, "ShearStiffness_H", ShearStiffness_H) 
    write(file, "NormalStiffness_H", NormalStiffness_H) 
    write(file, "Admissible", Admissible) 
    write(file, "Ranks_Shear", Ranks_Shear) 
    write(file, "Ranks_Normal", Ranks_Normal) 
    write(file, "ElementRange_SR", ElementRange_SR) 
else
    write(file, "StiffnessMatrixShear", StiffnessMatrixShear) 
    write(file, "StiffnessMatrixNormal", StiffnessMatrixNormal) 

end
close(file)

######################## Save Vertices for Triangles ##############################
if RorT == "T"        
    file = jldopen(OutputFileName, "a+")
    write(file, "P1", P1) 
    write(file, "P2", P2) 
    write(file, "P3", P3) 
    close(file)
end







