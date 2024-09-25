
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

function hmat_speed_test()

    GC.gc(false)
    SoluationHMatShear = 0.0
    SoluationOrignialShear = 0.0
    SoluationHMatNormal = 0.0
    SoluationOrignialNormal = 0.0

        TestVector = rand(FaultCount)

        for i=1:1
            Elastic_Load_EachThread = zeros(FaultCount, ThreadCount)
            SoluationTest = HmatSolver_Pararllel(TestVector, ShearStiffness_H, ElementRange_SR, 
                            Par_ElementDivision_Shear, ThreadCount, zeros(FaultCount, ThreadCount))
        # println("end")
        end

        ElapseTime_HMatShear = @elapsed for i=1:Repeats
            Elastic_Load_EachThread = zeros(FaultCount, ThreadCount)
            SoluationHMatShear = HmatSolver_Pararllel(TestVector, ShearStiffness_H, ElementRange_SR, 
                                                Par_ElementDivision_Shear, ThreadCount, zeros(FaultCount, ThreadCount))
            # println("end")
            end

        ElapseTime_HMatNormal = @elapsed for i=1:Repeats
            Elastic_Load_EachThread = zeros(FaultCount, ThreadCount)
            SoluationHMatNormal = HmatSolver_Pararllel(TestVector, NormalStiffness_H, ElementRange_SR, 
                                                Par_ElementDivision_Normal, ThreadCount, zeros(FaultCount, ThreadCount))
        # println("end")
        end


        if SaveOriginalMatrix == 1 
            ElapseTime_OriginalMatShear = @elapsed for i=1:Repeats
                SoluationOrignialShear = StiffnessMatrixShear * TestVector
            end

            ElapseTime_OrigianlMatNormal = @elapsed for i=1:Repeats
                SoluationOrignialNormal = StiffnessMatrixNormal * TestVector
            end
        end 
            # plot(SoluationHMatShear-SoluationOrignialShear)
            
        if SaveOriginalMatrix == 1 
            println("\nResults")
            println("Shear Product Hmat: ", ElapseTime_HMatShear, " Original Matrix: ", ElapseTime_OriginalMatShear, " Ratio: ", ElapseTime_HMatShear/ElapseTime_OriginalMatShear)
            println("Normal Product Hmat: ", ElapseTime_HMatNormal, " Original Matrix: ", ElapseTime_OrigianlMatNormal, " Ratio: ", ElapseTime_HMatNormal/ElapseTime_OrigianlMatNormal)
            println("Hmat Total: ", ElapseTime_HMatShear+ElapseTime_HMatNormal, " Original Matrix Total: ", ElapseTime_OriginalMatShear+ElapseTime_OrigianlMatNormal, 
                        " Ratio: ", (ElapseTime_HMatShear+ElapseTime_HMatNormal)/(ElapseTime_OriginalMatShear+ElapseTime_OrigianlMatNormal))
        else 
            println("\nResults")
            println("Shear Product Hmat: ", ElapseTime_HMatShear)
            println("Normal Product Hmat: ", ElapseTime_HMatNormal)
        end
    # plot(SoluationHMatShear - SoluationOrignialShear)
    
    GC.gc(true)
end

hmat_speed_test()

# sizeof(X)/1e6