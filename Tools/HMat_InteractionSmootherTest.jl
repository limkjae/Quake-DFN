
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
@pyimport matplotlib.patches as patches
pygui(true)

include("../Functions_Solvers.jl")
include("../Functions_RSFDFN3DMain_H.jl")
include("../Results/Functions_Plot.jl")
include("../QuickParameterAdjust.jl")
include("../Functions_Hmatrix.jl")
LoadingInputFileName="Input_Discretized.jld2" 


function ReduceInteractionforHmat()
    SaveIT = 1
    PlotIt = 1

    StrongInteractionCriteriaMultiple = 0.5


    # jldsave("HmatSave.jld2"; StiffnessMatrixShearOriginal, StiffnessMatrixNormalOriginal, Ranks, ElementRange_SR)
    Admissible = load(LoadingInputFileName,"Admissible")
    FaultCount= load(LoadingInputFileName, "FaultCount")
    Ranks_Shear= load(LoadingInputFileName, "Ranks_Shear") # figure(11); plot(Ranks)
    Ranks_Normal= load(LoadingInputFileName, "Ranks_Normal") # figure(11); plot(Ranks)
    ElementRange_SR = load(LoadingInputFileName, "ElementRange_SR")
    ShearStiffness_H = load(LoadingInputFileName, "ShearStiffness_H")
    NormalStiffness_H = load(LoadingInputFileName, "NormalStiffness_H")
    NormalStiffnessZero = load(LoadingInputFileName, "NormalStiffnessZero")
    # StiffnessMatrixShear= load(LoadingInputFileName, "StiffnessMatrixShear")
    # StiffnessMatrixNormal= load(LoadingInputFileName, "StiffnessMatrixNormal")


    LoadingStiffnessH, K_Self= StiffnessTransitionToLoading(ShearStiffness_H, ElementRange_SR, FaultCount)
    # typeof(LoadingStiffnessH[10]) == PartialQR{Float64}
    # typeof(LoadingStiffnessH[1877]) == PartialQR{Float64}

    # Admissible[Block] == 0
    # ElementRange_SR[Block,:]
    # K_Self
    # LoadingStiffnessH[Block]
    StrongInteractionPair = zeros(Int,1,2)
    # ElementRange_SR[Block,:]
    BlockCount = length(Admissible)
    for Block = 1:BlockCount
        if Admissible[Block] == 0
            SourceInThisBlock = 0
            for SourceAt = ElementRange_SR[Block,1]:ElementRange_SR[Block,2]
                SourceInThisBlock = SourceInThisBlock+1
                
                ReceiverInThisBlock = 0
                for ReceiverAt = ElementRange_SR[Block,3]:ElementRange_SR[Block,4]
                    ReceiverInThisBlock = ReceiverInThisBlock + 1
                    if abs(K_Self[ReceiverAt]) * StrongInteractionCriteriaMultiple < abs( LoadingStiffnessH[Block][ReceiverInThisBlock,SourceInThisBlock]) 
                    # println(SourceInThisBlock,"  ", ReceiverInThisBlock,"  ", K_Self[ReceiverAt], "  ", LoadingStiffnessH[Block][ReceiverInThisBlock,SourceInThisBlock])
                    ShearStiffness_H[Block][ReceiverInThisBlock,SourceInThisBlock] = 
                                sign(LoadingStiffnessH[Block][ReceiverInThisBlock,SourceInThisBlock]) * abs(K_Self[ReceiverAt]) * StrongInteractionCriteriaMultiple
                    StrongInteractionPair = [StrongInteractionPair;  [SourceAt ReceiverAt]]
                    end

                    # if 
                    # if LoadingStiffnessH[]
                
                end
            end
        end
    end
    StrongInteractionPair = StrongInteractionPair[2:end,:]

    println("Strong interaction Count: ", length(StrongInteractionPair[:,1]))


    if SaveIT == 1
        file = jldopen(LoadingInputFileName, "a+")

        Base.delete!(file, "ShearStiffness_H") 
        write(file, "ShearStiffness_H", ShearStiffness_H) 
        if length(StrongInteractionPair[:,1]) > 0
            println("ShearStiffness_H Adjusted")
        end
        close(file)
    end


    if PlotIt ==1
            
        FaultCount= load(LoadingInputFileName, "FaultCount")
        FaultCenter= load(LoadingInputFileName, "FaultCenter")
        FaultLengthStrike= load(LoadingInputFileName, "FaultLengthStrike")
        FaultLengthDip= load(LoadingInputFileName, "FaultLengthDip")
        FaultStrikeAngle= load(LoadingInputFileName, "FaultStrikeAngle")
        FaultDipAngle= load(LoadingInputFileName, "FaultDipAngle")
        FaultRakeAngle= load(LoadingInputFileName, "FaultRakeAngle")
        LoadingFaultCount= load(LoadingInputFileName, "LoadingFaultCount")
        Fault_a= load(LoadingInputFileName, "Fault_a")
        Fault_b= load(LoadingInputFileName, "Fault_b")
        Fault_Dc= load(LoadingInputFileName, "Fault_Dc")
        Fault_Theta_i= load(LoadingInputFileName, "Fault_Theta_i")
        Fault_V_i= load(LoadingInputFileName, "Fault_V_i")
        Fault_Friction_i= load(LoadingInputFileName, "Fault_Friction_i")
        Fault_NormalStress= load(LoadingInputFileName, "Fault_NormalStress")
        Fault_BulkIndex= load(LoadingInputFileName, "Fault_BulkIndex")
        # FaultMass= load(LoadingInputFileName, "FaultMass")


        InputProperty = ones(Int,FaultCount); ColorMinMax=[0,2]
        for i in eachindex(InputProperty)
            for j = 1:length(StrongInteractionPair[:,1])
                if i == StrongInteractionPair[j,1]
                    InputProperty[i] = 0
                elseif i == StrongInteractionPair[j,2]
                    InputProperty[i] = 2
                end
            end
        end

        #Single element locator
        # figure(1)
        SelectedElements = StrongInteractionPair[:,1]
        PlotRotation = [23,-16]
        Transparent = 1 # 1 for transparent fault plot. 0 for no-transparency
        Edge = 0 # 0 for no element boudary. 1 for plotting element boundary
        MinMax_Axis = 0 # 0 for automatically selected axis minimim and maximum 
        LoadingFaultPlot = 0 # 1 to plot constant velocity faults. 
        # clf()
        MaxVaule, MinValue = FaultPlot_3D_Color_AdjustedElements(FaultCenter,FaultLengthStrike, FaultLengthDip, FaultStrikeAngle,
                    FaultDipAngle, FaultRakeAngle, InputProperty, PlotRotation, MinMax_Axis, ColorMinMax,
                    Transparent, Edge, LoadingFaultCount)
    end
    
end

ReduceInteractionforHmat()