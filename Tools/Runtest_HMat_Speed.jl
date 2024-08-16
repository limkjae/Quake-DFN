
using PyPlot
using PyCall
using DelimitedFiles
using JLD2
using LinearAlgebra
using Printf
using SpecialFunctions
using HMatrices
using StaticArrays
using LowRankApprox
using Distributed
using LoopVectorization
using Statistics
pygui(true)

include("../Functions_Solvers.jl")
include("../Functions_RSFDFN3DMain_H.jl")
include("../Results/Functions_Plot.jl")
include("../QuickParameterAdjust.jl")
include("../Functions_Hmatrix.jl")
LoadingInputFileName="Input_Discretized.jld2" 

Repeats = 10

FaultCount= load(LoadingInputFileName, "FaultCount")
Ranks_Shear= load(LoadingInputFileName, "Ranks_Shear") # figure(11); plot(Ranks)
Ranks_Normal= load(LoadingInputFileName, "Ranks_Normal") # figure(11); plot(Ranks)
ElementRange_SR = load(LoadingInputFileName, "ElementRange_SR")
ShearStiffness_H = load(LoadingInputFileName, "ShearStiffness_H")
NormalStiffness_H = load(LoadingInputFileName, "NormalStiffness_H")
NormalStiffnessZero = load(LoadingInputFileName, "NormalStiffnessZero")
StiffnessMatrixShear= load(LoadingInputFileName, "StiffnessMatrixShear")
StiffnessMatrixNormal= load(LoadingInputFileName, "StiffnessMatrixNormal")




################ Optimizing and Initializing Hmatrix ################
# ThreadCount = Threads.nthreads()
ThreadCount = 16
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

function hmat_speed_test()



    TestVector = rand(FaultCount)

    SoluationHMatShear = 0.0
    SoluationOrignialShear = 0.0
    SoluationHMatNormal = 0.0
    SoluationOrignialNormal = 0.0

    Elastic_Load_EachThread = zeros(FaultCount, ThreadCount)
    for i=1:1
        SoluationTest = HmatSolver_Pararllel(TestVector, ShearStiffness_H, ElementRange_SR, 
                        Par_ElementDivision_Shear, ThreadCount, zeros(FaultCount, ThreadCount))
    # println("end")
    end

    HMatShear = @elapsed for i=1:Repeats
        Elastic_Load_EachThread = zeros(FaultCount, ThreadCount)
        SoluationHMatShear = HmatSolver_Pararllel(TestVector, ShearStiffness_H, ElementRange_SR, 
                                            Par_ElementDivision_Shear, ThreadCount, zeros(FaultCount, ThreadCount))
        # println("end")
        end

    OriginalMatShear = @elapsed for i=1:Repeats
        SoluationOrignialShear = StiffnessMatrixShear * TestVector
    end

    HMatNormal = @elapsed for i=1:Repeats
        Elastic_Load_EachThread = zeros(FaultCount, ThreadCount)
        SoluationHMatNormal = HmatSolver_Pararllel(TestVector, NormalStiffness_H, ElementRange_SR, 
                                            Par_ElementDivision_Normal, ThreadCount, zeros(FaultCount, ThreadCount))
    # println("end")
    end

    OrigianlMatNormal = @elapsed for i=1:Repeats
        SoluationOrignialNormal = StiffnessMatrixNormal * TestVector
    end
 
    println("\nResults")
    println("Shear Product Hmat: ", HMatShear, " Original Matrix: ", OriginalMatShear)
    println("Normal Product Hmat: ", HMatNormal, " Original Matrix: ", OrigianlMatNormal)
    # plot(SoluationHMatShear - SoluationOrignialShear)
end

hmat_speed_test()





# ################################################################## 
# NetDisp = rand(FaultCount)
# Elastic_Load_D = zeros(FaultCount)
# # Elastic_Load_EachThread = zeros(FaultCount, ThreadCount)
# Elastic_Load_D_SingleThread = zeros(FaultCount)
# Elastic_Load_D_FullMatrix  = zeros(FaultCount)
#     # Threads.@threads for ThreadIdx=1:ThreadCount 
# Elastic_Load_EachThread = zeros(FaultCount, ThreadCount)
        
# Repeats=10
# println("Parallel 1")




# for i=1:10
# println()
# @time Elastic_Load_D = HmatSolver_Pararllel(NetDisp, ShearStiffness_H, ElementRange_SR, 
#          Par_ElementDivision_Shear, ThreadCount, zeros(FaultCount, ThreadCount))
#         #  println(Elastic_Load_EachThread[1,1])
# # @time Elastic_Load_D = HmatSolver_Pararllel(NetDisp, ShearStiffness_H, ElementRange_SR, 
# # FaultCount, Par_ElementDivision_Shear, ThreadCount)

# @time Elastic_Load_D_FullMatrix = StiffnessMatrixShear * NetDisp

# plot(Elastic_Load_D - Elastic_Load_D_FullMatrix)
# end