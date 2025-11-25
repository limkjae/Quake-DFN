
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

println("Current Threads using: ",Threads.nthreads())

include("../../scripts/Functions_Solvers.jl")
include("../../scripts/Functions_RSFDFN3DMain_H.jl")
include("../../scripts/Functions_Plot.jl")
include("../../QuickParameterAdjust.jl")
include("../../scripts/Functions_Hmatrix.jl")
LoadingInputFileName="Input_Discretized.jld2" 


FaultCount= load(LoadingInputFileName, "FaultCount")
Ranks_Shear= load(LoadingInputFileName, "Ranks_Shear") # figure(11); plot(Ranks)
Admissible= load(LoadingInputFileName, "Admissible") # figure(11); plot(Ranks)
Ranks_Normal= load(LoadingInputFileName, "Ranks_Normal") # figure(11); plot(Ranks)
ElementRange_SR = load(LoadingInputFileName, "ElementRange_SR")
ShearStiffness_H = load(LoadingInputFileName, "ShearStiffness_H")
NormalStiffness_H = load(LoadingInputFileName, "NormalStiffness_H")
NormalStiffnessZero = load(LoadingInputFileName, "NormalStiffnessZero")


################ Optimizing and Initializing Hmatrix ################
Repeats = 100
ThreadCount = 8
Epsilon_MaxDiffRatio = 1e-7
MaxRatioAllowed = 1.5
MaxIteration = 50

BlockCount = length(Ranks_Shear)



# println("Shear Stiffness")
# Par_ElementDivision_Shear = ParallelOptimization(ShearStiffness_H, ElementRange_SR, 
#                             FaultCount, BlockCount, ThreadCount, MaxRatioAllowed, MaxIteration)
# println("Normal Stiffness")
# Par_ElementDivision_Normal = ParallelOptimization(NormalStiffness_H, ElementRange_SR, 
#                             FaultCount, BlockCount, ThreadCount, MaxRatioAllowed, MaxIteration)



println("Shear Stiffness")
Par_ElementDivision_Shear = ParallelOptimization2(ShearStiffness_H,  
                            BlockCount, ThreadCount, Ranks_Shear, Admissible)
Par_ElementDivision_Normal = Par_ElementDivision_Shear

# println("Normal Stiffness")
# Par_ElementDivision_Normal = ParallelOptimization2(NormalStiffness_H,  
#                             BlockCount, ThreadCount, Ranks_Normal)
# println("Normal Stiffness")
# Par_ElementDivision_Normal = ParallelOptimization2(NormalStiffness_H, ElementRange_SR, 
#                             FaultCount, BlockCount, ThreadCount, MaxRatioAllowed, MaxIteration)





# Complexity_S = zeros(BlockCount)
# Complexity_N = zeros(BlockCount)
# for BlockIdx=1:BlockCount
#     if Ranks_Shear[BlockIdx] > 0
#         Complexity_S[BlockIdx] = (size(ShearStiffness_H[BlockIdx])[1] + size(ShearStiffness_H[BlockIdx])[2]) * Ranks_Shear[BlockIdx]
#     else
#         Complexity_S[BlockIdx] = (size(ShearStiffness_H[BlockIdx])[1] * size(ShearStiffness_H[BlockIdx])[2]) 
#     end
    
#     if Ranks_Normal[BlockIdx] > 0
#         Complexity_N[BlockIdx] = (size(NormalStiffness_H[BlockIdx])[1] + size(NormalStiffness_H[BlockIdx])[2]) * Ranks_Normal[BlockIdx]
#     else
#         Complexity_N[BlockIdx] = (size(NormalStiffness_H[BlockIdx])[1] * size(NormalStiffness_H[BlockIdx])[2]) 
#     end
# end

# CumSum_S = cumsum(Complexity_S) 
# Division = CumSum_S[end] / ThreadCount
# Par_ElementDivision_Shear = [0]
# CurrentDiv = 1
# for BlockIdx=1:BlockCount

#     if CumSum_S[BlockIdx] > Division * CurrentDiv
#         Par_ElementDivision_Shear = [Par_ElementDivision_Shear; BlockIdx]
#         CurrentDiv += 1
#     end 

# end

# Par_ElementDivision_Shear = [Par_ElementDivision_Shear; BlockCount]
println(Par_ElementDivision_Shear)
println(Par_ElementDivision_Normal)

#####################################################################

function hmat_speed_test()

    GC.gc(false)
    SoluationHMatShear = 0.0
    SoluationOrignialShear = 0.0
    SoluationHMatNormal = 0.0
    SoluationOrignialNormal = 0.0

        # TestVector = rand(FaultCount)
        TestVector = ones(FaultCount)

        for i=1:1
            Elastic_Load_EachThread = zeros(FaultCount, ThreadCount)
            SoluationTest = HmatSolver_Pararllel(TestVector, ShearStiffness_H, ElementRange_SR, 
                            Par_ElementDivision_Shear, ThreadCount, zeros(FaultCount, ThreadCount))
        # println("end")
        end

        ElapseTime_HMatShear = @elapsed for i=1:Repeats
            
        # @time for i=1:Repeats
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


            # plot(SoluationHMatShear-SoluationOrignialShear)
            println(SoluationHMatShear[1:10])
            println("\nResults")
            println("Shear Product Hmat: ", ElapseTime_HMatShear)
            println("Normal Product Hmat: ", ElapseTime_HMatNormal)
    
    GC.gc(true)
end

hmat_speed_test()

# sizeof(X)/1e6