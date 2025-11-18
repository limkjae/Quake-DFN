
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

include("../../scripts/Functions_Solvers.jl")
include("../../scripts/Functions_RSFDFN3DMain_H.jl")
include("../../scripts/Functions_Plot.jl")
include("../../QuickParameterAdjust.jl")
include("../../scripts/Functions_Hmatrix.jl")
LoadingInputFileName="Input_Discretized.jld2" 

Repeats = 5

FaultCount= load(LoadingInputFileName, "FaultCount")
Ranks_Shear= load(LoadingInputFileName, "Ranks_Shear") # figure(11); plot(Ranks)
Ranks_Normal= load(LoadingInputFileName, "Ranks_Normal") # figure(11); plot(Ranks)
ElementRange_SR = load(LoadingInputFileName, "ElementRange_SR")
ShearStiffness_H = load(LoadingInputFileName, "ShearStiffness_H")
NormalStiffness_H = load(LoadingInputFileName, "NormalStiffness_H")
NormalStiffnessZero = load(LoadingInputFileName, "NormalStiffnessZero")




################ Optimizing and Initializing Hmatrix ################
ThreadREPL = Threads.nthreads()
ThreadCountAll = collect(1:1:8)
Epsilon_MaxDiffRatio = 1e-7
MaxRatioAllowed = 1.5
MaxIteration = 50

#####################################################################
# figure(1)
# clf()

function hmat_speed_test()
    
    SoluationOrignialShear = 0.0
    SoluationOrignialNormal = 0.0
    ElapseTime_HMatShear=zeros(length(ThreadCountAll), Repeats)
    ElapseTime_HMatNormal=zeros(length(ThreadCountAll), Repeats)
    ElapseTime_AllCalculation = zeros(length(ThreadCountAll), Repeats)

    GC.gc(false)
    TestVector = rand(FaultCount)

    GC.gc(true)

    for Tindex in eachindex(ThreadCountAll)
        ThreadCount = ThreadCountAll[Tindex]
        BlockCount = length(Ranks_Shear)
        println("Shear Stiffness")
        Par_ElementDivision_Shear = ParallelOptimization(ShearStiffness_H, ElementRange_SR, 
                                    FaultCount, BlockCount, ThreadCount, MaxRatioAllowed, MaxIteration)
        println("Normal Stiffness")
        Par_ElementDivision_Normal = ParallelOptimization(NormalStiffness_H, ElementRange_SR, 
                                    FaultCount, BlockCount, ThreadCount, MaxRatioAllowed, MaxIteration)
        
        SoluationHMatShear = 0.0
        SoluationHMatNormal = 0.0
        for Repindex = 1: Repeats

            GC.gc(false)
            for i=1:1
                Elastic_Load_EachThread = zeros(FaultCount, ThreadCount)
                SoluationTest = HmatSolver_Pararllel(TestVector, ShearStiffness_H, ElementRange_SR, 
                                Par_ElementDivision_Shear, ThreadCount, zeros(FaultCount, ThreadCount))
            # println("end")
            end

            ElapseTime_HMatShear[Tindex,Repindex] = @elapsed begin
                Elastic_Load_EachThread = zeros(FaultCount, ThreadCount)
                SoluationHMatShear = HmatSolver_Pararllel(TestVector, ShearStiffness_H, ElementRange_SR, 
                                                    Par_ElementDivision_Shear, ThreadCount, zeros(FaultCount, ThreadCount))
                # println("end")
                end


            ElapseTime_HMatNormal[Tindex,Repindex] = @elapsed begin
                Elastic_Load_EachThread = zeros(FaultCount, ThreadCount)
                SoluationHMatNormal = HmatSolver_Pararllel(TestVector, NormalStiffness_H, ElementRange_SR, 
                                                    Par_ElementDivision_Normal, ThreadCount, zeros(FaultCount, ThreadCount))
            # println("end")
            end
            GC.gc(true)
        end
                
        # println("\nResults")
        # println("Shear Product Hmat: ", ElapseTime_HMatShear, " Original Matrix: ", ElapseTime_OriginalMatShear, " Ratio: ", ElapseTime_HMatShear/ElapseTime_OriginalMatShear)
        # println("Normal Product Hmat: ", ElapseTime_HMatNormal, " Original Matrix: ", ElapseTime_OrigianlMatNormal, " Ratio: ", ElapseTime_HMatNormal/ElapseTime_OrigianlMatNormal)
        # println("Hmat Total: ", ElapseTime_HMatShear+ElapseTime_HMatNormal, " Original Matrix Total: ", ElapseTime_OriginalMatShear+ElapseTime_OrigianlMatNormal, 
        #             " Ratio: ", (ElapseTime_HMatShear+ElapseTime_HMatNormal)/(ElapseTime_OriginalMatShear+ElapseTime_OrigianlMatNormal))
        # # plot(SoluationHMatShear - SoluationOrignialShear)        
    end

    TotalElapsedTime = ElapseTime_HMatShear + ElapseTime_HMatNormal
    TotalAverage = sum(TotalElapsedTime, dims=2) / Repeats
    figure(1)
    clf()
    plot(ThreadCountAll, TotalElapsedTime,"ko",markerfacecolor="none" )
    plot(ThreadCountAll,TotalAverage,"k")
        
    MaxTimeInPlot = maximum(TotalAverage)
    plot([ThreadREPL,ThreadREPL],[0,MaxTimeInPlot*1.1],"r")
    ylim([0, MaxTimeInPlot*1.2])

end

hmat_speed_test()

