

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

function BuildInputFromBulkGeometry_H()

    OutputFileName="Input_Discretized.jld2"

    ##########################################################################
    ########################## Hmatrix compress? #############################
    HMatrixCompress = 1 # If this is 1, stiffness Matrix will be compressed by Hmatrix
    SaveOriginalMatrix = 1 # 1: save Original Matrix (can be very large), 0: Discard Original Matrix. 
    
    #####----- Hmatrix compression detail ----#####
    Tolerance = 1e3 # pascal for 1m slip (More approximaion for higher Tolerance)
  
    #####---------   HMatrix Plots?  --------#####
    PlotHMat = 1 # HMatrix structure plot
    PlotBlock3D = 0 # 3D block ploat
    ##########################################################################
    ##########################################################################
    

    ElementPartRoughCount = 2000
    
    @load "HmatrixStructure.jld2"
    FaultCount = length(Input_Segment[:,1])



    #####                               3D Plot                                   #####
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

    #####                     Hierarchical Matrix Structure Plot                  #####
    if PlotHMat == 1
        figure(2)
        println("Plotting Hmatrix Structure")
        clf()
        ax = gca()
        # ax[:set_aspect]("equal")
        for i=1:length(ElementRange_SR[:,1])
            if Admissible[i] > 0   
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
            


    ####################################################################################
    ############################# Build Stiffness Matrix ###############################

    StiffnessMatrixShearOriginal=zeros(FaultCount,FaultCount)
    StiffnessMatrixNormalOriginal=zeros(FaultCount,FaultCount)

    if Switch_StrikeSlip_or_ReverseNormal == 1
        println("Preparing for discretization")
        StiffnessMatrixShearOriginal, StiffnessMatrixNormalOriginal = 
            BuildMatrixByPartsShear(FaultCount, ElementPartRoughCount, Input_Segment,  ShearModulus, PoissonRatio)
        
    elseif Switch_StrikeSlip_or_ReverseNormal == 2
        println("Preparing for discretization")
        StiffnessMatrixShearOriginal, StiffnessMatrixNormalOriginal = 
            BuildMatrixByPartsNormal(FaultCount, ElementPartRoughCount, Input_Segment,  ShearModulus, PoissonRatio)

    else
        println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        println("!!!!!!! Switch_StrikeSlip_or_ReverseNormal should be either 1 or 2  !!!!!!!")
        println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    end

    # println(StiffnessMatrixNormalOriginal)

    NormalStiffnessZero = 0
    if iszero(StiffnessMatrixNormalOriginal)
        NormalStiffnessZero = 1
        println("Stiffness Matrix Normal is all zero")
    end
    

    ###################################################################################
    ############################### HMatrix Approximate ###############################
    if HMatrixCompress ==1
        BlockCount = length(ElementRange_SR[:,1])
        ShearStiffness_H = Any[0]
        Ranks_Shear = zeros(Int, BlockCount)
        NormalStiffness_H = Any[0]
        Ranks_Normal = zeros(Int, BlockCount)
        BlockIndex = 0
        println("compressing")
        for i=1:BlockCount
            BlockIndex = BlockIndex + 1
            
            if Admissible[BlockIndex] > 0
                OrigianlMatrixToApproximate = StiffnessMatrixShearOriginal[ElementRange_SR[i,1]:ElementRange_SR[i,2],ElementRange_SR[i,3]:ElementRange_SR[i,4]]
                ApproxMatrixS = pqrfact(OrigianlMatrixToApproximate, atol = Tolerance)
                push!(ShearStiffness_H,ApproxMatrixS)
                Ranks_Shear[BlockIndex] = size(ApproxMatrixS[:Q],2)
                
                OrigianlMatrixToApproximate = StiffnessMatrixNormalOriginal[ElementRange_SR[i,1]:ElementRange_SR[i,2],ElementRange_SR[i,3]:ElementRange_SR[i,4]]
                ApproxMatrixN = pqrfact(OrigianlMatrixToApproximate, atol = Tolerance)
                push!(NormalStiffness_H,ApproxMatrixN)
                Ranks_Normal[BlockIndex] = size(ApproxMatrixN[:Q],2)
            else 
                push!(ShearStiffness_H,StiffnessMatrixShearOriginal[ElementRange_SR[i,1]:ElementRange_SR[i,2],ElementRange_SR[i,3]:ElementRange_SR[i,4]])
                push!(NormalStiffness_H,StiffnessMatrixNormalOriginal[ElementRange_SR[i,1]:ElementRange_SR[i,2],ElementRange_SR[i,3]:ElementRange_SR[i,4]])
            end

        end
        ShearStiffness_H = ShearStiffness_H[2:end]
        NormalStiffness_H = NormalStiffness_H[2:end]

    end
    ########################## end of HMatrix Approximation ##########################


    # ########################## Remove Unstable Faults ##############################    
    # ReducedStiffnessMatrixShear, ReducedStiffnessMatrixNormal, ReducedInput_Segment=
    # CheckTooClose(StiffnessMatrixShearOriginal, StiffnessMatrixNormalOriginal, Input_Segment, Input_Bulk, DropCrit, DropCritNormalStressMultiplier);


    ################################### Save Files #################################
    # jldsave("HmatSave.jld2"; ShearStiffness_H, StiffnessMatrixShearOriginal, StiffnessMatrixNormalOriginal, Ranks, ElementRange_SR)
    if HMatrixCompress ==1
        SaveResults_H(StiffnessMatrixShearOriginal, StiffnessMatrixNormalOriginal, Input_Segment, NormalStiffnessZero,
            OutputFileName, ShearModulus, PoissonRatio, RockDensity, Switch_StrikeSlip_or_ReverseNormal, MinimumNS,
            Ranks_Shear, Ranks_Normal, ElementRange_SR, ShearStiffness_H, NormalStiffness_H, Admissible, SaveOriginalMatrix)

    else 
        SaveResults(StiffnessMatrixShearOriginal, StiffnessMatrixNormalOriginal, Input_Segment, NormalStiffnessZero,
            OutputFileName, ShearModulus, PoissonRatio, RockDensity, Switch_StrikeSlip_or_ReverseNormal, MinimumNS);
    end
end


Input = BuildInputFromBulkGeometry_H()
