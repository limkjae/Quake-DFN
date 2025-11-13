
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

include("../scripts/Functions_Solvers.jl")
include("../scripts/Functions_RSFDFN3DMain_H.jl")
include("../Results/Functions_Plot.jl")
include("../QuickParameterAdjust.jl")
include("../scripts/Functions_Hmatrix.jl")
LoadingInputFileName="Input_Discretized.jld2" 


PlotIt = 1
StrongInteractionCriteriaMultiple = 0.5
clf()

function CheckStrongInteractionHmat(StrongInteractionCriteriaMultiple, PlotIt)

    Admissible = load(LoadingInputFileName,"Admissible")
    FaultCount= load(LoadingInputFileName, "FaultCount")
    ElementRange_SR = load(LoadingInputFileName, "ElementRange_SR")
    ShearStiffness_H = load(LoadingInputFileName, "ShearStiffness_H")
    LoadingFaultCount= load(LoadingInputFileName, "LoadingFaultCount")

    LoadingStiffnessH, K_Self= StiffnessTransitionToLoading(ShearStiffness_H, ElementRange_SR, FaultCount)

    StrongInteractionPair = zeros(Int,1,2)

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
                        StrongInteractionPair = [StrongInteractionPair;  [SourceAt ReceiverAt]]
                    end

                end
            end
        end
    end
    StrongInteractionPair = StrongInteractionPair[2:end,:]
    println("HMatrix checked")
    println("Strong interaction Count: ", length(StrongInteractionPair[:,1]))


    if PlotIt ==1
        println("Plotting Element Pair")
        FaultCount= load(LoadingInputFileName, "FaultCount")
        FaultCenter= load(LoadingInputFileName, "FaultCenter")
        FaultLengthStrike= load(LoadingInputFileName, "FaultLengthStrike")
        FaultLengthDip= load(LoadingInputFileName, "FaultLengthDip")
        FaultStrikeAngle= load(LoadingInputFileName, "FaultStrikeAngle")
        FaultDipAngle= load(LoadingInputFileName, "FaultDipAngle")
        FaultRakeAngle= load(LoadingInputFileName, "FaultRakeAngle")
        LoadingFaultCount= load(LoadingInputFileName, "LoadingFaultCount")

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





function CheckStrongInteraction(StrongInteractionCriteriaMultiple, PlotIt)

    StiffnessMatrixShear= load(LoadingInputFileName, "StiffnessMatrixShear")
    Faultcount = load(LoadingInputFileName, "FaultCount")
    LoadingFaultCount= load(LoadingInputFileName, "LoadingFaultCount")
    #  SourceAt = 1
    #  ReceiverAt = 90
    StrongInteractionPair = zeros(Int,1,2)
    SourceInThisBlock = 0
    for SourceAt = 1:Faultcount - LoadingFaultCount
        SourceInThisBlock = SourceInThisBlock+1
        for ReceiverAt = 1:Faultcount - LoadingFaultCount
            if SourceAt != ReceiverAt
                if abs(StiffnessMatrixShear[ReceiverAt, ReceiverAt]) * StrongInteractionCriteriaMultiple < abs( StiffnessMatrixShear[ReceiverAt, SourceAt]) 
                    # println("Strong !!")
                    StrongInteractionPair = [StrongInteractionPair;  [SourceAt ReceiverAt]]
                end

            end
        end
    end
    StrongInteractionPair = StrongInteractionPair[2:end,:]
    println("Regular Matrix checked")
    println("Storng Interaction Pairs Count: ", length(StrongInteractionPair[:,1]))

    if PlotIt ==1
        println("Plotting Element Pair")
        FaultCount= load(LoadingInputFileName, "FaultCount")
        FaultCenter= load(LoadingInputFileName, "FaultCenter")
        FaultLengthStrike= load(LoadingInputFileName, "FaultLengthStrike")
        FaultLengthDip= load(LoadingInputFileName, "FaultLengthDip")
        FaultStrikeAngle= load(LoadingInputFileName, "FaultStrikeAngle")
        FaultDipAngle= load(LoadingInputFileName, "FaultDipAngle")
        FaultRakeAngle= load(LoadingInputFileName, "FaultRakeAngle")
        LoadingFaultCount= load(LoadingInputFileName, "LoadingFaultCount")

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





if haskey(load(LoadingInputFileName), "Admissible")
    CheckStrongInteractionHmat(StrongInteractionCriteriaMultiple, PlotIt)
else 
    CheckStrongInteraction(StrongInteractionCriteriaMultiple, PlotIt)
end





