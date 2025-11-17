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

############################################################################
############## Hmatrix with Distributed discretization #####################
Hmatrix = true
HowManyDistribution = 5
############################################################################


TriAlreadyCompiled = @isdefined FunctionRead
if TriAlreadyCompiled == false
    include("scripts/Functions_TDstressHS.jl")
    FunctionRead = 1
end

include("scripts/Functions_BuildInputFile.jl")
include("scripts/Functions_OKADA3D.jl")
include("scripts/Functions_Plot.jl")
include("scripts/Functions_Hmatrix.jl")




function DiscritizeFullMatrix()

    InputBulkFileName="Input_BulkFaultGeometry.txt"
    OutputFileName="Input_Discretized.jld2"
    HMatrixStructureFile = "Input_HmatrixStructure.jld2"



    Input_Segment, LoadingFaultCount, ShearModulus, PoissonRatio, RockDensity, 
    Switch_StrikeSlip_or_ReverseNormal, DropCrit, DropCritNormalStressMultiplier, MinimumNS, RorT, FaultCount =
        ReadBulkInput(InputBulkFileName)


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
        StiffnessMatrixShear, StiffnessMatrixNormal = 
            BuildMatrixByParts_Even_Rec(FaultCount, Input_Segment,  ShearModulus, PoissonRatio)
    elseif RorT == "T"
        StiffnessMatrixShear, StiffnessMatrixNormal = 
            BuildMatrixByParts_Even_Tri(P1, P2, P3, FaultRakeAngle, FaultCenter, UnitVector_Normal, UnitVector_Slip,
                                FaultCount, ShearModulus, PoissonRatio)                            
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
        write(file, "StiffnessMatrixShear", StiffnessMatrixShear) 
        write(file, "StiffnessMatrixNormal", StiffnessMatrixNormal) 
    close(file)

    ######################## Save Vertices for Triangles ##############################
    if RorT == "T"        
        file = jldopen(OutputFileName, "a+")
        write(file, "P1", P1) 
        write(file, "P2", P2) 
        write(file, "P3", P3) 
        close(file)
    end


end


if Hmatrix == true
    include("scripts/DistributedDiscretization.jl")

elseif Hmatrix == false
    DiscritizeFullMatrix()
end




