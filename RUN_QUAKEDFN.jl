

using PyPlot
using PyCall
using DelimitedFiles
using JLD2
using LinearAlgebra
using Printf
using SpecialFunctions
using StaticArrays
pygui(true)

include("Functions_Solvers.jl")
include("Functions_RSFDFN3DMain_D.jl")
include("Results/Functions_Plot.jl")
include("QuickParameterAdjust.jl")

SaveResultFileName="Results/Result.jld2"
SaveInputInfoFileName="Results/Result_Input.jld2" 

LoadingInputFileName="Input_Discretized.jld2" 


########################## Simulation Time Set ################################
TotalStep = 3000 # Total simulation step
SaveStep = 3000 # Automatically saved every this step
RecordStep = 10 # Simulation sampling rate


########################## Time Stepping Setup ################################
TimeStepOnlyBasedOnUnstablePatch = 1 # if 1, time step is calculated only based on the unstable patch
TimeStepPreset = 3 # 1: conservative --> 4: optimistic
RuptureTimeStepMultiple = 3
VerticalLengthScaleforM = 0 # if 0, Mass is automatically determined based on the fault length (radiation damping dominated). If not, M = VerticalLengthScaleforM * density / 2

# Manually adjust time step below. No change when 0.0
TimeSteppingAdj =   
        [0.0  0.0  0.0  0.0;   # Time step size
         0.0  0.0  0.0  0.0]   # Velocity

########## Strong Interaction Supression for Numerical Stability ##############
StrongInteractionCriteriaMultiple = 0.5 # only applied when larger than 0.

############################# Plots before run? ################################
DtPlot = 0 # 1 will plot dt vs maxV
GeometryPlot = 0 # 1 will plot a-b




function RunRSFDFN3D(TotalStep, RecordStep, 
    LoadingInputFileName, SaveResultFileName, RuptureTimeStepMultiple)


    ############################### Load Input Files ###############################
    ######++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++######
    StiffnessMatrixShear= load(LoadingInputFileName, "StiffnessMatrixShear")
    StiffnessMatrixNormal= load(LoadingInputFileName, "StiffnessMatrixNormal")
    FaultCenter= load(LoadingInputFileName, "FaultCenter")
    ShearModulus= load(LoadingInputFileName, "ShearModulus")
    RockDensity= load(LoadingInputFileName, "RockDensity")
    PoissonRatio= load(LoadingInputFileName, "PoissonRatio")
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
    MinimumNormalStress = load(LoadingInputFileName, "MinimumNormalStress")
    NormalStiffnessZero = load(LoadingInputFileName, "NormalStiffnessZero")
    ########^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^########
    ################################################################################


    if StrongInteractionCriteriaMultiple > 0    
        StiffnessMatrixShear, StiffnessMatrixNormal = 
            ReduceTooStrongInteraction(StrongInteractionCriteriaMultiple, FaultCount, StiffnessMatrixShear, StiffnessMatrixNormal)
    end



    if GeometryPlot==1
        PlotRotation=[45,-45]
        Transparent=0
        MinMax_Axis=[-3000 3000; -3000 3000; -4000 0]
        ColorMinMax=0
        PlotInput=Fault_a - Fault_b    
        Edge = 1
        clf()
        MaxVaule, MinValue = FaultPlot_3D_Color_General(FaultCenter,FaultLengthStrike, FaultLengthDip,
            FaultStrikeAngle, FaultDipAngle, FaultRakeAngle, PlotInput, 
            PlotRotation, MinMax_Axis, ColorMinMax, Transparent, Edge, LoadingFaultCount)
        figure(1).canvas.draw()
    end
    
    ############# Evolution and Alpha Value ################
    Alpha_Evo = 0.0
    EvolutionDR = 1
    ########################################################


    ################################ Run Simulation ################################
    ######++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++######


    ######+++++++++++++++++++++++++ Adjust Parameters ++++++++++++++++++++++++######
    if VerticalLengthScaleforM == 0
        VertScale = minimum([minimum(FaultLengthStrike), minimum(FaultLengthDip)])
    else 
        VertScale = VerticalLengthScaleforM
    end    
    FaultMass = ones(FaultCount) .* VertScale * RockDensity / 2


    LoadingFaultCount, FaultMass, Fault_a, Fault_b, Fault_Dc, Fault_Theta_i, Fault_V_i, 
    Fault_Friction_i, Fault_NormalStress, Fault_V_Const,  FaultCenter, FaultIndex_Adjusted, MinimumNormalStress = 
        ParameterAdj(LoadingFaultCount, FaultMass, Fault_a, Fault_b, Fault_Dc, Fault_Theta_i, Fault_V_i, 
        Fault_Friction_i, Fault_NormalStress, Fault_V_Const, 
        FaultStrikeAngle, FaultDipAngle, FaultCenter, Fault_BulkIndex, FaultRakeAngle, MinimumNormalStress)
    


    ######+++++++++++++++++++++++++   Time Step Set   ++++++++++++++++++++++++######
    Period=zeros(FaultCount)
    for i=1:FaultCount
        Period[i]=sqrt(FaultMass[i]/abs(StiffnessMatrixShear[i,i]))
    end
    RecTimeStep=minimum(Period)/10 
    println("Recommended TimeStep: ",RecTimeStep)
    RuptureDt = RecTimeStep* RuptureTimeStepMultiple

    if TimeStepPreset ==1
        global TimeStepping =
        [1e4 1e1 RuptureDt RuptureDt;
        1e-7 1e-5  1e-3 1e-2]

    elseif TimeStepPreset ==2
        global TimeStepping =
        [1e5 1e1 RuptureDt RuptureDt;
        1e-8 1e-5  1e-3 1e-2]
                
    elseif TimeStepPreset ==3
        global TimeStepping =
        [1e6 RecTimeStep*1000 RuptureDt RuptureDt;
        1e-9 1e-5  1e-3 1e-2]

    elseif TimeStepPreset ==4
        global TimeStepping =
        [1e7 RecTimeStep*1000 RuptureDt RuptureDt;
        1e-9 1e-5  1e-3 1e-2]

    end

    for i=1:length(TimeStepping[:,1])
        for j=1:length(TimeStepping[1,:])
            if TimeSteppingAdj[i,j] != 0
                global TimeStepping[i,j]=TimeSteppingAdj[i,j]
            end
        end
    end
    SwitchV=TimeStepping[2,3]
    if DtPlot==1
    FunctionDTPlot(SwitchV , TimeStepping, RecTimeStep)
    end

    if isfile(SaveResultFileName)
        println(" *******************************************************")
        println(" *** A Result File with Same Save Name Already Exist ***")
        println(" ***           This Can Generate Errors              ***")
        println(" *******************************************************")
    end

    ######+++++++++++++++++++++++++ Save Input Files ++++++++++++++++++++++++######
    save(SaveInputInfoFileName, 
    # "StiffnessMatrixShear", StiffnessMatrixShear, "StiffnessMatrixNormal", StiffnessMatrixNormal, 
    "FaultCenter", FaultCenter, "ShearModulus", ShearModulus, "RockDensity", RockDensity, "PoissonRatio", PoissonRatio,
    "FaultLengthStrike", FaultLengthStrike, "FaultLengthDip", FaultLengthDip, "FaultStrikeAngle", FaultStrikeAngle, 
    "FaultDipAngle", FaultDipAngle, "FaultRakeAngle", FaultRakeAngle, "Fault_a", Fault_a, "Fault_b", Fault_b, "Fault_Dc", Fault_Dc, 
    "Fault_Theta_i", Fault_Theta_i, "Fault_V_i", Fault_V_i, "Fault_Friction_i", Fault_Friction_i, "Fault_NormalStress", Fault_NormalStress, 
    "Fault_V_Const", Fault_V_Const, "Fault_BulkIndex", Fault_BulkIndex, "FaultLengthStrike_Bulk", FaultLengthStrike_Bulk, 
    "FaultLengthDip_Bulk", FaultLengthDip_Bulk, "FaultCount", FaultCount, "LoadingFaultCount", LoadingFaultCount, "FaultMass", FaultMass, "MinimumNormalStress", MinimumNormalStress)
    

    ######+++++++++++++++++++++++++         Run       ++++++++++++++++++++++++######
    main(StiffnessMatrixShear, StiffnessMatrixNormal, NormalStiffnessZero,
    ShearModulus, FaultCount, LoadingFaultCount, FaultMass,
    Fault_a, Fault_b, Fault_Dc, Fault_Theta_i, Fault_V_i, Fault_Friction_i,
    Fault_NormalStress, Fault_V_Const,
    TotalStep, RecordStep, SwitchV, TimeStepping, SaveResultFileName,RockDensity,
    FaultCenter,FaultLengthStrike, FaultLengthDip, FaultStrikeAngle, FaultDipAngle, FaultRakeAngle, SaveStep,
    TimeStepOnlyBasedOnUnstablePatch, MinimumNormalStress, Alpha_Evo, EvolutionDR)  


    ########^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^########
    ################################################################################

end


RunRSFDFN3D(TotalStep, RecordStep, 
    LoadingInputFileName, SaveResultFileName, RuptureTimeStepMultiple)

