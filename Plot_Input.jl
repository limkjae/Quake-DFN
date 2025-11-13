

using PyPlot
using PyCall
using DelimitedFiles
using JLD2
using LinearAlgebra
using Printf
using SpecialFunctions
using StaticArrays
using LowRankApprox
pygui(true)

include("scripts/Functions_Solvers.jl")
include("scripts/Functions_RSFDFN3DMain_D.jl")
include("Results/Functions_Plot.jl")
include("QuickParameterAdjust.jl")
include("scripts/Functions_Hmatrix.jl")

LoadingInputFileName="Input_Discretized.jld2" #only Needed when BuildStiffnessMatrix=2


function RunPlotInput(LoadingInputFileName)

    ################################################################################
    ############################### Load Input Files ###############################
    # StiffnessMatrixShear= load(LoadingInputFileName, "StiffnessMatrixShear")
    # StiffnessMatrixNormal= load(LoadingInputFileName, "StiffnessMatrixNormal")
    RorT= load(LoadingInputFileName, "RorT")
    FaultCenter= load(LoadingInputFileName, "FaultCenter")
    FaultLengthStrike= load(LoadingInputFileName, "FaultLengthStrike")
    FaultLengthDip= load(LoadingInputFileName, "FaultLengthDip")
    FaultStrikeAngle= load(LoadingInputFileName, "FaultStrikeAngle")
    FaultDipAngle= load(LoadingInputFileName, "FaultDipAngle")
    FaultRakeAngle= load(LoadingInputFileName, "FaultRakeAngle")
    Fault_a= load(LoadingInputFileName, "Fault_a")
    Fault_b= load(LoadingInputFileName, "Fault_b")
    Fault_Dc= load(LoadingInputFileName, "Fault_Dc")
    Fault_Theta_i= load(LoadingInputFileName, "Fault_Theta_i")
    Fault_V_i= load(LoadingInputFileName, "Fault_V_i")
    Fault_Friction_i= load(LoadingInputFileName, "Fault_Friction_i")
    Fault_NormalStress= load(LoadingInputFileName, "Fault_NormalStress")
    Fault_V_Const= load(LoadingInputFileName, "Fault_V_Const")
    Fault_BulkIndex= load(LoadingInputFileName, "Fault_BulkIndex")
    FaultLengthStrike_Bulk= load(LoadingInputFileName, "FaultLengthStrike_Bulk")
    FaultLengthDip_Bulk= load(LoadingInputFileName, "FaultLengthDip_Bulk")
    FaultCount= load(LoadingInputFileName, "FaultCount")
    LoadingFaultCount= load(LoadingInputFileName, "LoadingFaultCount")
    # FaultMass= load(LoadingInputFileName, "FaultMass")
    FaultMass = ones(FaultCount)
    MinimumNormalStress = load(LoadingInputFileName, "MinimumNormalStress")
    ########^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^########
    ################################################################################


    println("FaultCount: ", FaultCount)

    
    ################################################################################
    ######++++++++++++++++++++++ Apply Adjust Parameters +++++++++++++++++++++######
    LoadingFaultCount, FaultMass, Fault_a, Fault_b, Fault_Dc, Fault_Theta_i, Fault_V_i, Fault_Friction_i,
    Fault_NormalStress, Fault_V_Const, FaultCenter, FaultIndex_Adjusted = 
        ParameterAdj(LoadingFaultCount, FaultMass, Fault_a, Fault_b, Fault_Dc, Fault_Theta_i, Fault_V_i, 
        Fault_Friction_i, Fault_NormalStress, Fault_V_Const, 
        FaultStrikeAngle, FaultDipAngle, FaultCenter, Fault_BulkIndex, FaultRakeAngle, MinimumNormalStress)
        
    ########^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^########
    ################################################################################

    K_Self=[]
    if haskey(load(LoadingInputFileName), "ShearStiffness_H")
        ShearStiffness_H = load(LoadingInputFileName, "ShearStiffness_H")
        ElementRange_SR = load(LoadingInputFileName, "ElementRange_SR")
        LoadingStiffnessH, K_Self= StiffnessTransitionToLoading(ShearStiffness_H, ElementRange_SR, FaultCount)
        K_Self = -K_Self
    else 

        K_Self= diag(load(LoadingInputFileName, "StiffnessMatrixShear"))

    end


    ################################################################################
    ######++++++++++++++ Check which element is underresolved ++++++++++++++++######
        KoverKC = zeros(FaultCount)
        UnderResolved = zeros(FaultCount)
        for i=1:FaultCount
            
            KoverKC[i] =  -K_Self[i]/((Fault_b[i] - Fault_a[i])*Fault_NormalStress[i]/Fault_Dc[i]) 
            if KoverKC[i] < 0
                KoverKC[i] = 100
            end
            # KCoverK[i] =  ((Fault_b[i] - Fault_a[i])*Fault_NormalStress[i]/Fault_Dc[i]) /-StiffnessMatrixShear[i,i]
            if KoverKC[i] < 1 
            UnderResolved[i] = 1
            end
        end
    ########^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^########
    ################################################################################





    
    ################################################################################
    ####################### Which parameter want to plot? ##########################
        ColorMinMax = 0  
        # PlotInput = log10.(Fault_Theta_i); ColorMinMax = 0 
        # PlotInput = log10.(Fault_V_i); ColorMinMax = 0  
        PlotInput =Fault_NormalStress; ColorMinMax = 0    
        # PlotInput =KoverKC ;ColorMinMax=[0,5]
        # PlotInput =UnderResolved ;ColorMinMax=[0,1]
        # PlotInput = Fault_a - Fault_b; ColorMinMax = 0  
        # PlotInput =  Fault_BulkIndex; ColorMinMax = 0  
        # PlotInput = Fault_Dc; ColorMinMax = 0  
        # PlotInput = FaultRakeAngle; ColorMinMax = 0  
        # PlotInput = FaultRakeAngle .* Fault_Friction_i; ColorMinMax = 0  
        # PlotInput = Fault_Friction_i; ColorMinMax = 0  

        PlotRotation=[45,-30]
        Edge = 0
        Transparent = 0
        MinMax_Axis=0 # automatically detect max and min 
        # MinMax_Axis=[-2000 2000; -2000 2000; -5000 0]
        # LoadingFaultCount=0
    ########^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^########
    ################################################################################
        





    ################################################################################
    ####################################### Plot ###################################
    if RorT == "R"
        figure(1)
        clf()
        MaxVaule, MinValue = FaultPlot_3D_Color_General(FaultCenter,FaultLengthStrike, FaultLengthDip,
            FaultStrikeAngle, FaultDipAngle, FaultRakeAngle, PlotInput, 
            PlotRotation, MinMax_Axis, ColorMinMax, Transparent, Edge, LoadingFaultCount)
        
        # figure(1)
        plotforcbar=  scatter([1,1],[1,1],0.1, [MinValue,MaxVaule], cmap="jet")
        colorbar(plotforcbar, pad=0.15)
        figure(1).canvas.draw()
        xlabel("x")
        ylabel("y")
    else
        
        if ColorMinMax == 0 
        MaxValue=maximum(PlotInput)
        MinValue=minimum(PlotInput)
        else
        MaxValue=ColorMinMax[2]
        MinValue=ColorMinMax[1]
        end

        figure(1)
        clf()
        art3d = PyObject(PyPlot.art3D)
        ax = subplot(projection="3d")
        if Edge == 0 
            edge_color = [0.2, 0.2, 0.2, 0.0]
        else
            edge_color = [0.2, 0.2, 0.2, 0.2]
        end

        for ElemIdx = 1:FaultCount- LoadingFaultCount
            cm = get_cmap(:jet)
            PlotValue=(PlotInput[ElemIdx]-MinValue)/(MaxValue-MinValue)

            if Transparent ==0
                face_color = [cm(PlotValue)[1], cm(PlotValue)[2],cm(PlotValue)[3],1.0]
            else
                face_color = [cm(PlotValue)[1], cm(PlotValue)[2],cm(PlotValue)[3],0.5]
            end

            verts = ((P1[ElemIdx,:],P2[ElemIdx,:],P3[ElemIdx,:]), )
            p3c = PyObject(art3d.Poly3DCollection(verts))
            pycall(ax.add_collection3d, PyAny, p3c)
            pycall(p3c.set_facecolor, PyAny, face_color)
            pycall(p3c.set_edgecolor, PyAny, edge_color)
            ax.view_init(PlotRotation[1], PlotRotation[2])
        end

        plotforcbar=  scatter([1,1],[1,1],0.1, [MinValue,MaxValue], cmap="jet")
        colorbar(plotforcbar, pad=0.15)
        figure(1).canvas.draw()
        ax.set_aspect("equal")

    end
    ########^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^########
    ################################################################################





    ################################################################################
    ############################ Single element locator ############################
    #=
    figure(1)
     # SelectedElements = FaultIndex_Adjusted
    SelectedElements = [843]
    PlotInput = Fault_a - Fault_b; ColorMinMax = 0  
        PlotRotation=[45,-30]
        Edge = 0
        Transparent = 0
        MinMax_Axis=0 # automatically detect max and min 
    clf()
     MaxVaule, MinValue = FaultPlot_3D_Color_SelectedElements(FaultCenter,FaultLengthStrike, FaultLengthDip,
         FaultStrikeAngle, FaultDipAngle, FaultRakeAngle, PlotInput, 
         PlotRotation, MinMax_Axis, ColorMinMax, Transparent, SelectedElements)
    =#
    ########^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^########
    ################################################################################
        
end


RunPlotInput(LoadingInputFileName)




