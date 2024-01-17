

using PyPlot
using PyCall
using DelimitedFiles
using JLD2
using LinearAlgebra
using Printf
using SpecialFunctions
using HMatrices
using StaticArrays
pygui(true)

include("Functions_Solvers.jl")
include("Functions_RSFDFN3DMain_D.jl")
include("Results/Functions_Plot.jl")
include("QuickParameterAdjust.jl")

LoadingInputFileName="Input_Discretized.jld2" #only Needed when BuildStiffnessMatrix=2


function RunPlotInput(LoadingInputFileName)

    ################################################################################
    ############################### Load Input Files ###############################
    StiffnessMatrixShear= load(LoadingInputFileName, "StiffnessMatrixShear")
    StiffnessMatrixNormal= load(LoadingInputFileName, "StiffnessMatrixNormal")
    FaultCenter= load(LoadingInputFileName, "FaultCenter")
    FaultLengthStrike= load(LoadingInputFileName, "FaultLengthStrike")
    FaultLengthDip= load(LoadingInputFileName, "FaultLengthDip")
    FaultStrikeAngle= load(LoadingInputFileName, "FaultStrikeAngle")
    FaultDipAngle= load(LoadingInputFileName, "FaultDipAngle")
    FaultLLRR= load(LoadingInputFileName, "FaultLLRR")
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
    FaultMass= load(LoadingInputFileName, "FaultMass")
    MinimumNormalStress = load(LoadingInputFileName, "MinimumNormalStress")
    ########^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^########
    ################################################################################




    
    ################################################################################
    ######++++++++++++++++++++++ Apply Adjust Parameters +++++++++++++++++++++######
    LoadingFaultCount, FaultMass, Fault_a, Fault_b, Fault_Dc, Fault_Theta_i, Fault_V_i, Fault_Friction_i,
    Fault_NormalStress, Fault_V_Const, StiffnessMatrixShear, StiffnessMatrixNormal, FaultCenter, FaultIndex_Adjusted = 
        ParameterAdj(LoadingFaultCount, FaultMass, Fault_a, Fault_b, Fault_Dc, Fault_Theta_i, Fault_V_i, 
        Fault_Friction_i, Fault_NormalStress, Fault_V_Const, StiffnessMatrixShear, StiffnessMatrixNormal, 
        FaultStrikeAngle, FaultDipAngle, FaultCenter, Fault_BulkIndex, FaultLLRR, MinimumNormalStress)
        
    ########^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^########
    ################################################################################

    



    ################################################################################
    ######++++++++++++++ Check which element is underresolved ++++++++++++++++######
        KoverKC = zeros(FaultCount)
        UnderResolved = zeros(FaultCount)
        for i=1:FaultCount
            
            KoverKC[i] =  -StiffnessMatrixShear[i,i]/((Fault_b[i] - Fault_a[i])*Fault_NormalStress[i]/Fault_Dc[i]) 
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
        # PlotInput =Fault_NormalStress; ColorMinMax = 0    
        # PlotInput =KoverKC ;ColorMinMax=[0,5]
        # PlotInput =UnderResolved ;ColorMinMax=[0,1]
        PlotInput = Fault_a - Fault_b; ColorMinMax = 0  
        # PlotInput =  Fault_BulkIndex; ColorMinMax = 0  
        # PlotInput = Fault_Dc; ColorMinMax = 0  
        # PlotInput = FaultLLRR; ColorMinMax = 0  
        # PlotInput = FaultLLRR .* Fault_Friction_i; ColorMinMax = 0  
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
    figure(1)
    clf()
    MaxVaule, MinValue = FaultPlot_3D_Color_General(FaultCenter,FaultLengthStrike, FaultLengthDip,
        FaultStrikeAngle, FaultDipAngle, FaultLLRR, PlotInput, 
        PlotRotation, MinMax_Axis, ColorMinMax, Transparent, Edge, LoadingFaultCount)
    
    # figure(1)
    plotforcbar=  scatter([1,1],[1,1],0.1, [MinValue,MaxVaule], cmap="jet")
    colorbar(plotforcbar, pad=0.15)
    figure(1).canvas.draw()
    xlabel("x")
    ylabel("y")
    ########^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^########
    ################################################################################





    ################################################################################
    ############################ Single element locator ############################
    
    # figure(1)
    # # SelectedElements = FaultIndex_Adjusted
    # SelectedElements = [20]
    # # clf()
    # MaxVaule, MinValue = FaultPlot_3D_Color_SelectedElements(FaultCenter,FaultLengthStrike, FaultLengthDip,
    #     FaultStrikeAngle, FaultDipAngle, FaultLLRR, PlotInput, 
    #     PlotRotation, MinMax_Axis, ColorMinMax, Transparent, SelectedElements)

    ########^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^########
    ################################################################################
        
end


RunPlotInput(LoadingInputFileName)







    
# ResultTime=load("Results/Result.jld","History_Time")
# ResultDisp=load("Results/Result.jld","History_Disp")
# ResultV=load("Results/Result.jld","History_V")

# figure(3)
# clf()
# PyPlot.plot(ResultTime[1:end-1], log10.(ResultV[1:end-1,:]), linewidth=1)
# figure(4)
# PyPlot.plot( log10.(ResultV[1:end-1,:]), linewidth=1)
