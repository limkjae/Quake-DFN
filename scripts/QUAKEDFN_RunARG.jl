
using PyPlot
using PyCall
using DelimitedFiles
using JLD2
using LinearAlgebra
using Printf
using SpecialFunctions
# using StaticArrays
using LowRankApprox
using Distributed
using Statistics
pygui(true)
GC.gc()

include("Functions_Solvers.jl")
include("Functions_Plot.jl")
include("../QuickParameterAdjust.jl")
include("Functions_Hmatrix.jl")
include("Functions_RSFDFN3DMain_H.jl")
include("Functions_RSFDFN3DMain_D.jl")



            if Sys.isapple()
                SaveResultFileName="$(@__DIR__)/../Results/$(ARGS[10]).jld2"
                SaveInputInfoFileName="$(@__DIR__)/../Results/$(ARGS[10])_Input.jld2" 
                LoadingInputFileName="$(@__DIR__)/../Input_Discretized.jld2" 
            else
                SaveResultFileName="Results/$(ARGS[10]).jld2"
                SaveInputInfoFileName="Results/$(ARGS[10])_Input.jld2" 
                LoadingInputFileName="Input_Discretized.jld2" 
            end


function RunRSFDFN3D()

    ########################## Simulation Time Set ################################
    TotalStep = parse(Int, ARGS[1]) # Total simulation step
    SaveStep = parse(Int, ARGS[2]) # Automatically saved every this step
    RecordStep = parse(Int, ARGS[3]) # Simulation sampling rate !! should be a factor of SaveStep !!

    ThreadCount = parse(Int, ARGS[4]) # if zero, it uses current thread count opened in REPL (Only for Hmatrix)

    ########################## Time Stepping Setup ################################
    DtCut = parse(Float64, ARGS[5])
    SwitchV = parse(Float64, ARGS[6])
    RuptureTimeStepMultiple =  parse(Float64, ARGS[7])
    MaximumDt =  parse(Float64, ARGS[8])
    VerticalLengthScaleforM =  parse(Float64, ARGS[9]) # if 0, Mass is automatically determined based on the fault length (radiation damping dominated for large rupture). If not, M = VerticalLengthScaleforM * density / 2


    ########## Strong Interaction Supression for Numerical Stability ##############
    StrongInteractionCriteriaMultiple = parse(Float64, ARGS[11]) # only applied when larger than 0. The higher, the more tolerance of strong interaction. 


    ############################# Plots before run? ################################
    DtPlot = 0 # 1 will plot dt vs maxV
    GeometryPlot = 0 # 1 will plot a-b

    Hmatrix = load(LoadingInputFileName, "Hmatrix")

    if Hmatrix == true
        

        println(" --- H Matrix Simulation --- ")
        if Threads.nthreads() == 1
            println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            println("          Only 1 thread is being used             ")
            println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        end
    else
        println("Full Matrix will be used")
    end

    ################################################################################
    ############################### Load Input Files ###############################

    RorT = load(LoadingInputFileName, "RorT")
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
        
    if Hmatrix == true
        Ranks_Shear= load(LoadingInputFileName, "Ranks_Shear") # figure(11); plot(Ranks)
        Ranks_Normal= load(LoadingInputFileName, "Ranks_Normal") # figure(11); plot(Ranks)
        ElementRange_SR = load(LoadingInputFileName, "ElementRange_SR")
        ShearStiffness_H = load(LoadingInputFileName, "ShearStiffness_H")
        NormalStiffness_H = load(LoadingInputFileName, "NormalStiffness_H")
        Admissible = load(LoadingInputFileName,"Admissible")
    else
        StiffnessMatrixShear= load(LoadingInputFileName, "StiffnessMatrixShear")
        StiffnessMatrixNormal= load(LoadingInputFileName, "StiffnessMatrixNormal")

    end

    if RorT == "T"
        P1 = load(LoadingInputFileName, "P1")
        P2 = load(LoadingInputFileName, "P2")
        P3 = load(LoadingInputFileName, "P3")
    end
    ################################################################################
    
    if Hmatrix == true
        if StrongInteractionCriteriaMultiple > 0
            ShearStiffness_H, NormalStiffness_H = ReduceTooStrongInteraction_Hmat(StrongInteractionCriteriaMultiple, Admissible,
            FaultCount, ElementRange_SR, ShearStiffness_H, NormalStiffness_H)
        end
    end
    # NormalStiffnessZero = 1

    ################################################################################
    ################################# Plot Geometry ################################
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
    ################################################################################



    

    ################################################################################
    ############################### Adjust Parameters ##############################

    #####----- Mass Calucation from element Size -----#####
    if VerticalLengthScaleforM == 0
        VertScale = minimum([minimum(FaultLengthStrike), minimum(FaultLengthDip)])
    else 
        VertScale = VerticalLengthScaleforM
    end    
    FaultMass = ones(FaultCount) .* VertScale * RockDensity / 2
    # FaultMass = minimum([FaultLengthStrike FaultLengthDip],dims = 2) * RockDensity ./ 3
    # FaultMass .= 1e6  


    ############# Evolution and Alpha Value ################
    Alpha_Evo = 0.0
    EvolutionDR = 1
    ########################################################
    
    #######---------- Logical quick adjust ----------#######
    LoadingFaultCount, FaultMass, Fault_a, Fault_b, Fault_Dc, Fault_Theta_i, Fault_V_i, 
    Fault_Friction_i, Fault_NormalStress, Fault_V_Const, FaultCenter, FaultIndex_Adjusted, MinimumNormalStress = 
        ParameterAdj(LoadingFaultCount, FaultMass, Fault_a, Fault_b, Fault_Dc, Fault_Theta_i, Fault_V_i, 
        Fault_Friction_i, Fault_NormalStress, Fault_V_Const, 
        FaultStrikeAngle, FaultDipAngle, FaultCenter, Fault_BulkIndex, FaultRakeAngle, MinimumNormalStress)

    ############################### Adjust Parameters ##############################
    ################################################################################
    


    ################################################################################
    ################################ Time Step Set  ################################

    Period=zeros(FaultCount)
    if Hmatrix == true
        K_Self  = GetKself(ShearStiffness_H, ElementRange_SR, FaultCount)
    else 
        K_Self = diag(StiffnessMatrixShear)
    end

    for i=1:FaultCount
        Period[i]=sqrt(FaultMass[i]/abs(K_Self[i]))
    end
    RecTimeStep=minimum(Period)/10 
    println("Recommended TimeStep: ",RecTimeStep)
    RuptureDt = RecTimeStep* RuptureTimeStepMultiple


    ################################ Time Step Set  ################################
    ################################################################################






    ###############################################################################
    ############################### Save Input Files ##############################
    save(SaveInputInfoFileName, 
    # "StiffnessMatrixShear", StiffnessMatrixShear, "StiffnessMatrixNormal", StiffnessMatrixNormal, 
    "FaultCenter", FaultCenter, "ShearModulus", ShearModulus, "RockDensity", RockDensity, "PoissonRatio", PoissonRatio,
    "FaultLengthStrike", FaultLengthStrike, "FaultLengthDip", FaultLengthDip, "FaultStrikeAngle", FaultStrikeAngle, 
    "FaultDipAngle", FaultDipAngle, "FaultRakeAngle", FaultRakeAngle, "Fault_a", Fault_a, "Fault_b", Fault_b, "Fault_Dc", Fault_Dc, 
    "Fault_Theta_i", Fault_Theta_i, "Fault_V_i", Fault_V_i, "Fault_Friction_i", Fault_Friction_i, "Fault_NormalStress", Fault_NormalStress, 
    "Fault_V_Const", Fault_V_Const, "Fault_BulkIndex", Fault_BulkIndex, "FaultLengthStrike_Bulk", FaultLengthStrike_Bulk, 
    "FaultLengthDip_Bulk", FaultLengthDip_Bulk, "FaultCount", FaultCount, "LoadingFaultCount", LoadingFaultCount, "FaultMass", FaultMass, 
    "MinimumNormalStress", MinimumNormalStress,"RorT",RorT, "Hmatrix", Hmatrix)
    
    if Hmatrix == true
        file = jldopen(SaveInputInfoFileName, "a+")
        # write(file, "ShearStiffness_H", ShearStiffness_H) 
        # write(file, "NormalStiffness_H", NormalStiffness_H) 
        write(file, "Admissible", Admissible) 
        write(file, "Ranks_Shear", Ranks_Shear) 
        write(file, "Ranks_Normal", Ranks_Normal) 
        write(file, "ElementRange_SR", ElementRange_SR) 
        close(file)
    end
            
    if RorT == "T"        
        file = jldopen(SaveInputInfoFileName, "a+")
        write(file, "P1", P1) 
        write(file, "P2", P2) 
        write(file, "P3", P3) 
        close(file)

    end
    

    ###############################################################################
    ################################ Run Simulation ###############################
    if Hmatrix == true
        main_H(ShearModulus, FaultCount, LoadingFaultCount, FaultMass, NormalStiffnessZero,
        Fault_a, Fault_b, Fault_Dc, Fault_Theta_i, Fault_V_i, Fault_Friction_i,
        Fault_NormalStress, Fault_V_Const,
        TotalStep, RecordStep, SwitchV, DtCut, RuptureDt, MaximumDt, SaveResultFileName,RockDensity,
        FaultCenter,FaultLengthStrike, FaultLengthDip, FaultStrikeAngle, FaultDipAngle, FaultRakeAngle, SaveStep,
        MinimumNormalStress, Alpha_Evo,
        Ranks_Shear, Ranks_Normal, ElementRange_SR, NormalStiffness_H, ShearStiffness_H, ThreadCount, EvolutionDR)    
    else
        main(StiffnessMatrixShear, StiffnessMatrixNormal, NormalStiffnessZero,
        ShearModulus, FaultCount, LoadingFaultCount, FaultMass,
        Fault_a, Fault_b, Fault_Dc, Fault_Theta_i, Fault_V_i, Fault_Friction_i,
        Fault_NormalStress, Fault_V_Const,
        TotalStep, RecordStep, SwitchV, DtCut, RuptureDt, MaximumDt, SaveResultFileName,RockDensity,
        FaultCenter,FaultLengthStrike, FaultLengthDip, FaultStrikeAngle, FaultDipAngle, FaultRakeAngle, SaveStep,
        MinimumNormalStress, Alpha_Evo, EvolutionDR)   
    end
    ################################ Run Simulation ###############################
    ###############################################################################



end


RunRSFDFN3D() 
println("Simulation finished. Press any key to close this window")
readline()


