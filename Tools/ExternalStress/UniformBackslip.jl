
using DelimitedFiles
using JLD2
using LinearAlgebra
using Printf
using SpecialFunctions
using PyPlot
using PyCall
using LowRankApprox
using Statistics
pygui(true)


LoadingRate=1e-9

OutputFile="Input_ExternalStressChange.jld2"
LoadingInputFileName="Input_Discretized.jld2" 

include("../../scripts/Functions_Hmatrix.jl")

############################### Load Input Files ###############################
######++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++######
RorT = load(LoadingInputFileName, "RorT")
FaultCenter= load(LoadingInputFileName, "FaultCenter")
ShearModulus= load(LoadingInputFileName, "ShearModulus")
RockDensity= load(LoadingInputFileName, "RockDensity")
PoissonRatio= load(LoadingInputFileName, "PoissonRatio")
FaultLengthStrike= load(LoadingInputFileName, "FaultLengthStrike")
FaultLengthDip= load(LoadingInputFileName, "FaultLengthDip")
FaultStrikeAngle= load(LoadingInputFileName, "FaultStrikeAngle")
FaultDipAngle= load(LoadingInputFileName, "FaultDipAngle")
FaultRakeAngle= load(LoadingInputFileName, "FaultRakeAngle")
Hmatrix = load(LoadingInputFileName, "Hmatrix")
FaultCount= load(LoadingInputFileName, "FaultCount")

if Hmatrix == true
    Ranks_Shear= load(LoadingInputFileName, "Ranks_Shear") # figure(11); plot(Ranks)
    Ranks_Normal= load(LoadingInputFileName, "Ranks_Normal") # figure(11); plot(Ranks)
    ElementRange_SR = load(LoadingInputFileName, "ElementRange_SR")
    ShearStiffness_H = load(LoadingInputFileName, "ShearStiffness_H")
    NormalStiffness_H = load(LoadingInputFileName, "NormalStiffness_H")
    Admissible = load(LoadingInputFileName,"Admissible")
else
    StiffnessMatrixShear= load(LoadingInputFileName, "StiffnessMatrixShear")
    StiffnessMatrixNormal= load(LoadingInputFileName, "StiffnessMatrixNormal")

end

if RorT == "T"
    P1 = load(LoadingInputFileName, "P1")
    P2 = load(LoadingInputFileName, "P2")
    P3 = load(LoadingInputFileName, "P3")
end


################ Optimizing and Initializing Hmatrix ################
Repeats = 10
ThreadCount = 5
Epsilon_MaxDiffRatio = 1e-7
MaxRatioAllowed = 1.5
MaxIteration = 50

BlockCount = length(Ranks_Shear)

println("Shear Stiffness")
Par_ElementDivision_Shear = ParallelOptimization(ShearStiffness_H, ElementRange_SR, 
                            FaultCount, BlockCount, ThreadCount, MaxRatioAllowed, MaxIteration)
println("Normal Stiffness")
Par_ElementDivision_Normal = ParallelOptimization(NormalStiffness_H, ElementRange_SR, 
                            FaultCount, BlockCount, ThreadCount, MaxRatioAllowed, MaxIteration)

########^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^########
################################################################################

StreeVectorShear = HmatSolver_Pararllel( ones(FaultCount), ShearStiffness_H, ElementRange_SR, 
                Par_ElementDivision_Shear, ThreadCount, zeros(FaultCount, ThreadCount))

StreeVectorNormal = HmatSolver_Pararllel( ones(FaultCount), NormalStiffness_H, ElementRange_SR, 
                Par_ElementDivision_Normal, ThreadCount, zeros(FaultCount, ThreadCount))

ExternalStress_Normal = zeros(2,FaultCount)
ExternalStress_Shear =  zeros(2,FaultCount)
ExternalStress_Shear[2,:] = - StreeVectorShear* 1e11 * 1e-9
ExternalStress_Normal[2,:] =  - StreeVectorNormal* 1e11 * 1e-9
ExternalStress_TimeArray = zeros(2)
ExternalStress_TimeArray[2] = 1e11
Pressure =  zeros(2,FaultCount)
FlowRate = 0
PressureOrigin = [0,0,0]


save(OutputFile, 
"ExternalStress_Normal", ExternalStress_Normal, "ExternalStress_Shear", ExternalStress_Shear, "Pressure", Pressure,
"ExternalStress_TimeArray", ExternalStress_TimeArray, "FlowRate", FlowRate, "PressureOrigin", PressureOrigin)
    
