
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

function hmat_speed_test()



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
    ThreadCount = Threads.nthreads()
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

    TestVector = rand(FaultCount)


    @elapsed for i=1:1
    SoluationHMat = HmatSolver_Pararllel(TestVector, ShearStiffness_H, ElementRange_SR, FaultCount,
                                        Par_ElementDivision_Shear, ThreadCount)
    # println("end")
    end

    HMatShear = @elapsed for i=1:Repeats
        SoluationHMat = HmatSolver_Pararllel(TestVector, ShearStiffness_H, ElementRange_SR, FaultCount,
                                            Par_ElementDivision_Shear, ThreadCount)
        # println("end")
        end
    OriginalMatShear = @elapsed for i=1:Repeats
    SoluationOrignial = StiffnessMatrixShear * TestVector
    end

    HMatNormal = @elapsed for i=1:Repeats
        SoluationHMat = HmatSolver_Pararllel(TestVector, NormalStiffness_H, ElementRange_SR, FaultCount,
                                            Par_ElementDivision_Shear, ThreadCount)
    # println("end")
    end

    OrigianlMatNormal = @elapsed for i=1:Repeats
        SoluationOrignial = StiffnessMatrixNormal * TestVector
    end
        
    # for i=1:5
    #     @time SoluationHMat = HmatSolver_ThreadTime(TestVector, ShearStiffness_H, ElementRange_SR, FaultCount,
    #                                          Par_ElementDivision_Shear, ThreadCount)
    #     println("done")
    # end

    println("\n Results")
    println("Shear Product Hmat: ", HMatShear, " Original Matrix: ", OriginalMatShear)
    println("Normal Product Hmat: ", HMatNormal, " Original Matrix: ", OrigianlMatNormal)

end

hmat_speed_test()