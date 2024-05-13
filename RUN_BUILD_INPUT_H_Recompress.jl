

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

include("Functions_BuildInputFile.jl")
include("Functions_OKADA3D.jl")
include("Results/Functions_Plot.jl")
include("Functions_Hmatrix.jl")

function Recompress()

    LoadingInputFileName="Input_Discretized.jld2" 
    OutputFileName="Input_Discretized.jld2"

    Tolerance = 1e3 # pascal for 1m slip

    ############# Plots? ##############
    PlotHMat = 1
    PlotBlock3D = 0
    ###################################


    StiffnessMatrixShearOriginal= load(LoadingInputFileName, "StiffnessMatrixShear")
    StiffnessMatrixNormalOriginal= load(LoadingInputFileName, "StiffnessMatrixNormal")
    Ranks= load(LoadingInputFileName, "Ranks") # figure(11); plot(Ranks)
    ElementRange_SR = load(LoadingInputFileName, "ElementRange_SR")
    ShearStiffness_H = load(LoadingInputFileName, "ShearStiffness_H")
    FaultCount =  load(LoadingInputFileName, "FaultCount")
    Admissible =  load(LoadingInputFileName, "Admissible")

    figure(10)
    plot(Ranks,"b")

    ################################################################################
    #################################### 3D Plot ###################################

    if PlotBlock3D == 1
        println("Plotting 3D Group")
        PlotRotation=[45,-30]
        Edge = 0
        Transparent = 0
        MinMax_Axis=0 # automatically detect max and min 
        LoadingFaultCount=0 
        ColorMinMax = 0  
        figure(1)
        clf()
        MaxVaule, MinValue = FaultPlot_3D_Color_General(Input_Segment[:,1:3],Input_Segment[:,4], Input_Segment[:,5],
            Input_Segment[:,6], Input_Segment[:,7], Input_Segment[:,8], Input_Segment[:,20], 
            PlotRotation, MinMax_Axis, ColorMinMax, Transparent, Edge, LoadingFaultCount)

        plotforcbar=  scatter([1,1],[1,1],0.1, [MinValue,MaxVaule], cmap="jet")
        colorbar(plotforcbar, pad=0.15)
        figure(1).canvas.draw()
        xlabel("x")
        ylabel("y")
    end
    ########^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^########
    ################################################################################



    ################################################################################
    ########################### Plot Hierarchical Matrix ###########################
    if PlotHMat == 1
        figure(2)
        println("Plotting Hmatrix Structure")
        clf()
        ax = gca()
        # ax[:set_aspect]("equal")
        for i=1:length(ElementRange_SR[:,1])
            if Ranks[i] > 0   
                c = PyObject(patches.Rectangle((ElementRange_SR[i,3]-1, -ElementRange_SR[i,1]+1), ElementRange_SR[i,4] - ElementRange_SR[i,3]+1, 
                            -ElementRange_SR[i,2] + ElementRange_SR[i,1]-1, linewidth=1, edgecolor="k", facecolor=[0.4  0.4  1]))    
            else           
                c = PyObject(patches.Rectangle((ElementRange_SR[i,3]-1, -ElementRange_SR[i,1]+1), ElementRange_SR[i,4] - ElementRange_SR[i,3]+1, 
                            -ElementRange_SR[i,2] + ElementRange_SR[i,1]-1, linewidth=1, edgecolor="k", facecolor=[1 0.4 0.4]))
            end
            # ax.text( (ElementRange_SR[i,3] + ElementRange_SR[i,4])/2, -(ElementRange_SR[i,1] + ElementRange_SR[i,2])/2, i ,size=8, horizontalalignment="center", verticalalignment="center", color="k") 
            ax.add_patch(c) 
        end
        xlim(0,FaultCount)
        ylim(-FaultCount,0)
    end
    ########^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^########
    ################################################################################



    ################################## Approximate ######################################
    BlockCount = length(ElementRange_SR[:,1])
    ShearStiffness_H = Any[0]
    BlockIndex = 0
    # Ranks = Ranks .* 2
    println("compressing")
    for i=1:BlockCount
        BlockIndex = BlockIndex + 1
        
        if Admissible[BlockIndex] > 0
            OrigianlMatrixToApproximate = StiffnessMatrixShearOriginal[ElementRange_SR[i,1]:ElementRange_SR[i,2],ElementRange_SR[i,3]:ElementRange_SR[i,4]]
            ApproxMatrix = pqrfact(OrigianlMatrixToApproximate, atol = Tolerance)
            # ApproxMatrix = pqrfact(OrigianlMatrixToApproximate, rank = 1)
            push!(ShearStiffness_H,ApproxMatrix)
            Ranks[BlockIndex] = size(ApproxMatrix[:Q],2)
            
        else 
            push!(ShearStiffness_H,StiffnessMatrixShearOriginal[ElementRange_SR[i,1]:ElementRange_SR[i,2],ElementRange_SR[i,3]:ElementRange_SR[i,4]])
        end

    end
    ShearStiffness_H = ShearStiffness_H[2:end]
    figure(10)
    plot(Ranks,"r")

    # ########################## Remove Unstable Faults ##############################    
    # ReducedStiffnessMatrixShear, ReducedStiffnessMatrixNormal, ReducedInput_Segment=
    # CheckTooClose(StiffnessMatrixShearOriginal, StiffnessMatrixNormalOriginal, Input_Segment, Input_Bulk, DropCrit, DropCritNormalStressMultiplier);


    ################################### Save Files #################################
    # jldsave("HmatSave.jld2"; ShearStiffness_H, StiffnessMatrixShearOriginal, StiffnessMatrixNormalOriginal, Ranks, ElementRange_SR)


    # @load "Input_Discretized.jld2"
    jldopen(OutputFileName, "a+") do file
        Base.delete!(file, "Ranks")
        Base.delete!(file, "ShearStiffness_H")
        write(file, "Ranks", Ranks)
        write(file, "ShearStiffness_H", ShearStiffness_H)
    end

    # SaveResults_H(StiffnessMatrixShearOriginal, StiffnessMatrixNormalOriginal, Input_Segment,
    #      OutputFileName, ShearModulus, PoissonRatio, RockDensity, Switch_StrikeSlip_or_ReverseNormal, MinimumNS,
    #      Ranks, ElementRange_SR, ShearStiffness_H)

end


Input = Recompress()
