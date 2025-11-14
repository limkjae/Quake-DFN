

function ParameterAdj(LoadingFaultCount, FaultMass, Fault_a, Fault_b, Fault_Dc, 
    Fault_Theta_i, Fault_V_i, Fault_Friction_i, Fault_NormalStress, Fault_V_Const,
     FaultStrikeAngle, FaultDipAngle, FaultCenter, Fault_BulkIndex, FaultRakeAngle, MinimumNormalStress)

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




    

    ######################################################################################################
    ##########################  Calculation of initial state from stress orientation #####################
    
    # MaxStressOrientation = 85. # between 0-180 degree
    # StressRatioMaxOverMin = 0.5
    # MinFrictionAllowed = 0.1 # smaller than this friction is not allowed

    # StressGradAtMaxOrientation = 6000.0
    # SurfaceStressAtMaxOrientation = 2e6
    # Fault_Theta_i .= 1e10
    # Fault_V_i .= 0.0
    # Friction_0 = ones(FaultCount) * 0.30
    # V0=1e-9;

    # Fault_Friction_i, Fault_NormalStress, Fault_V_i, Fault_Theta_i = 
    #             StressDependentFrictionParametersStrikeSlip(MaxStressOrientation, StressRatioMaxOverMin, MinFrictionAllowed,
    #             StressGradAtMaxOrientation, SurfaceStressAtMaxOrientation,
    #             FaultStrikeAngle, FaultDipAngle, Fault_V_i, Fault_Theta_i, Fault_Friction_i, 
    #             Fault_a, Fault_b, Fault_Dc, Fault_NormalStress, Friction_0, FaultCenter)

    ##########^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#############
    ######################################################################################################
    
    

    ######################################################################################################
    ######################################### Direct Adjust ##############################################
    # for i in eachindex(Fault_Dc)
    #     if Fault_BulkIndex[i]==1
    #         Fault_V_i[i]=1e-12
    #     end
    # end
    
    # for i=1:FaultCount
    #     if -1500 > FaultCenter[i,3] || FaultCenter[i,3] > -100
    #         Fault_a[i] = 0.01
    #     end
    # end

    # Fault_Theta_i .= 0.9e10
    # Fault_Dc .= 2e-3
    # Fault_a .= 0.05
    # Fault_b .= 0.003
    # Fault_NormalStress .= 10e6
    # Fault_V_i .= Fault_V_i * 10
    # Fault_Theta_i .= 1e10
    # MinimumNormalStress= 1e6

    # FaultMass .= 1e5
    # FaultMass .= 1e6
    ###^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^###
    #####################################################################################################


    if FaultMass_Original != FaultMass
        println("- Fault Mass Adjusted")
    end
    
    if Fault_a_Original != Fault_a
        println("- Fault a Adjusted")
    end
    if Fault_b_Original != Fault_b
        println("- Fault b Adjusted")
    end
    if Fault_Dc_Original != Fault_Dc
        println("- Fault Dc Adjusted")
    end
    if Fault_Theta_i_Original != Fault_Theta_i
        println("- Fault Theta_i Adjusted")
    end
    if Fault_V_i_Original != Fault_V_i
        println("- Fault Fault_V_i Adjusted")
    end
    if Fault_Friction_i_Original != Fault_Friction_i
        println("- Fault Friction_i Adjusted")
    end
    if Fault_NormalStress_Original != Fault_NormalStress
        println("- Fault Normal stress Adjusted")
    end
    if Fault_V_Const_Original != Fault_V_Const
        println("- Fault V_Const Adjusted")
    end
    if FaultCenter_Original != FaultCenter
        println("- Fault Center Adjusted")
    end

    if MinimumNormalStress_Original != MinimumNormalStress
        println("- Minimum NormalStress Adjusted")
    end

    return LoadingFaultCount, FaultMass, Fault_a, Fault_b, Fault_Dc, Fault_Theta_i, Fault_V_i, 
    Fault_Friction_i, Fault_NormalStress, Fault_V_Const, FaultCenter, FaultIndex_Adjusted, MinimumNormalStress


end



function StressDependentFrictionParametersStrikeSlip(MaxStressOrientation, StressRatioMaxOverMin, MinFrictionAllowed,
    StressGradAtMaxOrientation, SurfaceStressAtMaxOrientation,
    FaultStrikeAngle, FaultDipAngle, Fault_V_i, Fault_Theta_i, Fault_Friction_i, 
    Fault_a, Fault_b, Fault_Dc, Fault_NormalStress, Friction_0, FaultCenter)
    V0 = 1e-9
    NormalStressParameter = (1+StressRatioMaxOverMin)/2 .+ (1-StressRatioMaxOverMin)/2 .* cosd.(2 .* (FaultStrikeAngle .- 90.0 .- MaxStressOrientation))
    ShearStressParameter = -(1-StressRatioMaxOverMin)/2 .* sind.(2 .* (FaultStrikeAngle .- 90.0 .- MaxStressOrientation))
    Fault_Friction_i .= abs.(ShearStressParameter ./ NormalStressParameter)
    
        
    for i in eachindex(Fault_Friction_i)
        if Fault_Friction_i[i] < MinFrictionAllowed
            Fault_Friction_i[i] = MinFrictionAllowed
        end
        Fault_NormalStress[i] = (StressGradAtMaxOrientation * FaultCenter[i,3] + SurfaceStressAtMaxOrientation) * NormalStressParameter[i] 
    end
    
    if iszero(Fault_V_i)
        Fault_V_i = V0 .* exp.( (Fault_Friction_i .- Friction_0 .- Fault_b .* log.(Fault_Theta_i .* V0./Fault_Dc)) ./ Fault_a)
    end
    
    
    if iszero(Fault_Theta_i)
        Fault_Theta_i = Fault_Dc ./ V0 .* exp.( (Fault_Friction_i .- Friction_0 .- Fault_a .* log.(Fault_V_i ./ V0)) ./ Fault_b)
    end
    
    return Fault_Friction_i, Fault_NormalStress, Fault_V_i, Fault_Theta_i
end
