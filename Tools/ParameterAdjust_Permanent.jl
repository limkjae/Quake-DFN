
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



LoadingInputFileName="Input_Discretized.jld2" 

FaultCenter= load(LoadingInputFileName, "FaultCenter")
FaultMass= load(LoadingInputFileName, "FaultMass")
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
LoadingFaultCount= load(LoadingInputFileName, "LoadingFaultCount")
MinimumNormalStress = load(LoadingInputFileName, "MinimumNormalStress")



function ParameterAdj_permanent(LoadingFaultCount, FaultMass, Fault_a, Fault_b, Fault_Dc, 
    Fault_Theta_i, Fault_V_i, Fault_Friction_i, Fault_NormalStress, Fault_V_Const,
     FaultStrikeAngle, FaultDipAngle, FaultCenter, Fault_BulkIndex, FaultLLRR, MinimumNormalStress)

    FaultMass_Original = copy(FaultMass)
    Fault_a_Original = copy(Fault_a)
    Fault_b_Original = copy(Fault_b)
    Fault_Dc_Original = copy(Fault_Dc)
    Fault_Theta_i_Original = copy(Fault_Theta_i)
    Fault_V_i_Original = copy(Fault_V_i)
    Fault_Friction_i_Original = copy(Fault_Friction_i)
    Fault_NormalStress_Original = copy(Fault_NormalStress)
    Fault_V_Const_Original = copy(Fault_V_Const)
    FaultCenter_Original = copy(FaultCenter)
    FaultCount = length(Fault_a)
    MinimumNormalStress_Original = copy(MinimumNormalStress)  
    Count=0; 
    FaultIndex_Adjusted=0




    if FaultIndex_Adjusted == 0
        println("Adjusted Stiffness Count: 0")
    else
        println("Adjusted Stiffness Count: ", length(FaultIndex_Adjusted))
    end

    ##########^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#############
    ######################################################################################################


    #####################################################################################################
    #########################  Calculation of initial state from stress orientation #####################
    
    # MaxStressOrientation = 150. # between 0-180 degree
    # StressRatioMaxOverMin = 0.5
    # MinFrictionAllowed = 0.1 # smaller than this friction is not allowed

    # StressGradAtMaxOrientation = 6000.0
    # SurfaceStressAtMaxOrientation = 2e6
    # Fault_Theta_i .= 1e10
    # Fault_V_i .= 0.0
    # Friction_0 = ones(FaultCount) * 0.32
    # V0=1e-9;

    # Fault_Friction_i, Fault_NormalStress, Fault_V_i, Fault_Theta_i = 
    #             StressDependentFrictionParameters(MaxStressOrientation, StressRatioMaxOverMin, MinFrictionAllowed,
    #             StressGradAtMaxOrientation, SurfaceStressAtMaxOrientation,
    #             FaultStrikeAngle, FaultDipAngle, Fault_V_i, Fault_Theta_i, Fault_Friction_i, FaultLLRR,
    #             Fault_a, Fault_b, Fault_Dc, Fault_NormalStress, Friction_0, FaultCenter)

    #########^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#############
    #####################################################################################################
    
    

    ######################################################################################################
    ######################################### Direct Adjust ##############################################
    # for i in eachindex(Fault_Dc)
    #     if Fault_BulkIndex[i]==1
    #         Fault_a[i]=0.003
    #     end
    # end
    
    # for i=1:FaultCount
    #     if  FaultCenter[i,3] < 500  ||  3000 < FaultCenter[i,3] 
    #         Fault_a[i] = 0.002
    #     end
    # end

    # Fault_Theta_i .= 1.5027018579405773e9
    # Fault_Dc .= 3e-3
    # Fault_a .= 0.05
    # Fault_b .= 0.003
    # Fault_NormalStress .= 10e6
    # Fault_V_i .= 1e-13
    # Fault_Theta_i .= 1e10
    # MinimumNormalStress= 1e6

    # FaultMass .= 1e6
    ###^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^###
    #####################################################################################################




    file = jldopen(LoadingInputFileName, "a+")

    if FaultMass_Original != FaultMass
        Base.delete!(file, "FaultMass") 
        write(file, "FaultMass", FaultMass) 
        println("- Fault Mass Adjusted")
    end
    
    if Fault_a_Original != Fault_a
        Base.delete!(file, "Fault_a") 
        write(file, "Fault_a", Fault_a) 
        println("- Fault a Adjusted")
    end
    if Fault_b_Original != Fault_b
        Base.delete!(file, "Fault_b") 
        write(file, "Fault_b", Fault_b) 
        println("- Fault b Adjusted")
    end
    if Fault_Dc_Original != Fault_Dc
        Base.delete!(file, "Fault_Dc") 
        write(file, "Fault_Dc", Fault_Dc) 
        println("- Fault Dc Adjusted")
    end
    if Fault_Theta_i_Original != Fault_Theta_i
        Base.delete!(file, "Fault_Theta_i") 
        write(file, "Fault_Theta_i", Fault_Theta_i) 
        println("- Fault Theta_i Adjusted")
    end
    if Fault_V_i_Original != Fault_V_i
        Base.delete!(file, "Fault_V_i") 
        write(file, "Fault_V_i", Fault_V_i) 
        println("- Fault Fault_V_i Adjusted")
    end
    if Fault_Friction_i_Original != Fault_Friction_i
        Base.delete!(file, "Fault_Friction_i") 
        write(file, "Fault_Friction_i", Fault_Friction_i) 
        println("- Fault Friction_i Adjusted")
    end
    if Fault_NormalStress_Original != Fault_NormalStress
        Base.delete!(file, "Fault_NormalStress") 
        write(file, "Fault_NormalStress", Fault_NormalStress) 
        println("- Fault Normal stress Adjusted")
    end
    if Fault_V_Const_Original != Fault_V_Const
        Base.delete!(file, "Fault_V_Const") 
        write(file, "Fault_V_Const", Fault_V_Const) 
        println("- Fault V_Const Adjusted")
    end
    if FaultCenter_Original != FaultCenter
        Base.delete!(file, "FaultCenter") 
        write(file, "FaultCenter", FaultCenter) 
        println("- Fault Center Adjusted")
    end 

    if MinimumNormalStress_Original != MinimumNormalStress
        Base.delete!(file, "MinimumNormalStress") 
        write(file, "MinimumNormalStress", MinimumNormalStress) 
        println("- Minimum NormalStress Adjusted")
    end

    close(file)

    return LoadingFaultCount, FaultMass, Fault_a, Fault_b, Fault_Dc, Fault_Theta_i, Fault_V_i, 
    Fault_Friction_i, Fault_NormalStress, Fault_V_Const, FaultCenter, FaultIndex_Adjusted, MinimumNormalStress


end



function StressDependentFrictionParameters(MaxStressOrientation, StressRatioMaxOverMin, MinFrictionAllowed,
    StressGradAtMaxOrientation, SurfaceStressAtMaxOrientation,
    FaultStrikeAngle, FaultDipAngle, Fault_V_i, Fault_Theta_i, Fault_Friction_i, FaultLLRR,
    Fault_a, Fault_b, Fault_Dc, Fault_NormalStress, Friction_0, FaultCenter)
    V0 = 1e-9
    NormalStressParameter = (1+StressRatioMaxOverMin)/2 .+ (1-StressRatioMaxOverMin)/2 .* cosd.(2 .* (FaultStrikeAngle .- 90.0 .- MaxStressOrientation))
    ShearStressParameter = -(1-StressRatioMaxOverMin)/2 .* sind.(2 .* (FaultStrikeAngle .- 90.0 .- MaxStressOrientation))
    Fault_Friction_i .= abs.(ShearStressParameter ./ NormalStressParameter)
    
    if FaultLLRR != sign.(ShearStressParameter)
        println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        println("!!! Warning: Stress Orientation inconsistant with Slip Direction !!! ")
        println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    end
        
    for i in eachindex(Fault_Friction_i)
        if Fault_Friction_i[i] < MinFrictionAllowed
            Fault_Friction_i[i] = MinFrictionAllowed
        end
        Fault_NormalStress[i] = StressGradAtMaxOrientation * NormalStressParameter[i] * FaultCenter[i,3] + SurfaceStressAtMaxOrientation
    end
    
    if iszero(Fault_V_i)
        Fault_V_i = V0 .* exp.( (Fault_Friction_i .- Friction_0 .- Fault_b .* log.(Fault_Theta_i .* V0./Fault_Dc)) ./ Fault_a)
    end
    
    
    if iszero(Fault_Theta_i)
        Fault_Theta_i = Fault_Dc ./ V0 .* exp.( (Fault_Friction_i .- Friction_0 .- Fault_a .* log.(Fault_V_i ./ V0)) ./ Fault_b)
    end
    
    return Fault_Friction_i, Fault_NormalStress, Fault_V_i, Fault_Theta_i
end



ParameterAdj_permanent(LoadingFaultCount, FaultMass, Fault_a, Fault_b, Fault_Dc, 
    Fault_Theta_i, Fault_V_i, Fault_Friction_i, Fault_NormalStress, Fault_V_Const,
     FaultStrikeAngle, FaultDipAngle, FaultCenter, Fault_BulkIndex, FaultLLRR, MinimumNormalStress);