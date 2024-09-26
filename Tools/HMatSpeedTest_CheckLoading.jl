
using PyPlot
using PyCall
using DelimitedFiles
using JLD2
using LinearAlgebra
using Printf
using SpecialFunctions
using StaticArrays
using LowRankApprox
using Distributed
using Statistics
pygui(true)

include("../Functions_Solvers.jl")
include("../Functions_RSFDFN3DMain_H.jl")
include("../Results/Functions_Plot.jl")
include("../QuickParameterAdjust.jl")
include("../Functions_Hmatrix.jl")
LoadingInputFileName="Input_Discretized.jld2" 


FaultCount= load(LoadingInputFileName, "FaultCount")
Ranks_Shear= load(LoadingInputFileName, "Ranks_Shear") # figure(11); plot(Ranks)
Ranks_Normal= load(LoadingInputFileName, "Ranks_Normal") # figure(11); plot(Ranks)
ElementRange_SR = load(LoadingInputFileName, "ElementRange_SR")
ShearStiffness_H = load(LoadingInputFileName, "ShearStiffness_H")
NormalStiffness_H = load(LoadingInputFileName, "NormalStiffness_H")
NormalStiffnessZero = load(LoadingInputFileName, "NormalStiffnessZero")
SaveOriginalMatrix =  load(LoadingInputFileName, "SaveOriginalMatrix")

if SaveOriginalMatrix == 1 
    StiffnessMatrixShear= load(LoadingInputFileName, "StiffnessMatrixShear")
    StiffnessMatrixNormal= load(LoadingInputFileName, "StiffnessMatrixNormal")
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

#####################################################################


    GC.gc(false)
    SoluationHMatShear = 0.0
    SoluationOrignialShear = 0.0
    SoluationHMatNormal = 0.0
    SoluationOrignialNormal = 0.0
        
    TestVector = rand(FaultCount)

    SoluationHMatShear = HmatSolver_Pararllel(TestVector, ShearStiffness_H, ElementRange_SR, 
                                        Par_ElementDivision_Shear, ThreadCount, zeros(FaultCount, ThreadCount))

        SoluationOrignialShear = StiffnessMatrixShear * TestVector

        Elastic_Load_EachThread = zeros(FaultCount, ThreadCount)
        @inbounds Threads.@threads for ThreadIdx in 1:ThreadCount 

            for Blockidx = Par_ElementDivision_Shear[ThreadIdx]+1:Par_ElementDivision_Shear[ThreadIdx+1]
                @inbounds @views Elastic_Load_EachThread[ElementRange_SR[Blockidx,3]:ElementRange_SR[Blockidx,4], ThreadIdx] +=
                ShearStiffness_H[Blockidx] * TestVector[ElementRange_SR[Blockidx,1]:ElementRange_SR[Blockidx,2]]
            end   
      
                
        end
        Elastic_Load_DispP = sum(Elastic_Load_EachThread, dims=2)

        plot(SoluationHMatShear - SoluationOrignialShear)
