function CalculateFrictionStress_R(BulkFaultCount, Input_Bulk, MinFrictionAllowed,StressRatioXYZ,LoadingFaultInvert,
    StressOnSurface_Sig1Orientation, StressGredient_Sig1Orientation, MinimumNormalStressAllowed)
        ShearFric = zeros(BulkFaultCount)
        NormalFric = zeros(BulkFaultCount)
    for BulkIndex = 1:BulkFaultCount

        FaultStrikeAngle = Input_Bulk[BulkIndex, 6]
        FaultDipAngle = Input_Bulk[BulkIndex, 7]
        # Rotate the fault center stress to reference frame to read shear and normal
        RotationMat_FromFault_Strike=
        [cosd(-FaultStrikeAngle) -sind(-FaultStrikeAngle)  0
        sind(-FaultStrikeAngle) cosd(-FaultStrikeAngle) 0
        0  0  1];
        RotationMat_FromFault_Dip=
        [1 0 0
        0 cosd(-FaultDipAngle) -sind(-FaultDipAngle)
        0 sind(-FaultDipAngle) cosd(-FaultDipAngle)]
        RotationMat_FromFault_All = RotationMat_FromFault_Dip * RotationMat_FromFault_Strike
        
        StressOnFault = RotationMat_FromFault_All * StressRatioXYZ * RotationMat_FromFault_All'    
        # RakeAngle[BulkIndex] = rad2deg(atan(StressOnFault[2,3] / StressOnFault[1,3]))
        Rake = rad2deg(atan(StressOnFault[2,3] / StressOnFault[1,3]))
        ShearFric[BulkIndex] = sqrt(StressOnFault[2,3]^2 + StressOnFault[1,3]^2) 
        NormalFric[BulkIndex] = abs(StressOnFault[3,3])
        
        Input_Bulk[BulkIndex,15] = NormalFric[BulkIndex] .* StressOnSurface_Sig1Orientation
        if Input_Bulk[BulkIndex,15] <  MinimumNormalStressAllowed
            Input_Bulk[BulkIndex,15] = MinimumNormalStressAllowed
        end
        Input_Bulk[BulkIndex,16] = NormalFric[BulkIndex] .* StressGredient_Sig1Orientation 
        Friction = sqrt(StressOnFault[2,3]^2 + StressOnFault[1,3]^2) / abs(StressOnFault[3,3])

        if Friction < MinFrictionAllowed
            println("Fault ", BulkIndex, " has very small friction. Adjusted to Minimum Allowed Friction")
            ShearFric[BulkIndex] = ShearFric[BulkIndex]  * MinFrictionAllowed / Friction
            Friction = MinFrictionAllowed
        end
        if Input_Bulk[BulkIndex,17] > 0 && LoadingFaultInvert == 1
            Rake = Rake + 180
        elseif Input_Bulk[BulkIndex,17] > 0 && LoadingFaultAdjust == 0
            # println("do nothing")
        else
            if StressOnFault[1,3] < 0;  Rake = Rake + 180;  end
            if Rake > 360; Rake = Rake - 360; end
            if Rake < 0; Rake = Rake + 360; end
            Input_Bulk[BulkIndex, 8] =  Rake
            Input_Bulk[BulkIndex, 14] =  Friction
        end
        # println(Rake, "   ", Friction)
    end
    return Input_Bulk
end

function CalculateFrictionStress_T(FaultCount, Input_Bulk, MinFrictionAllowed,StressRatioXYZ , 
    StressOnSurface_Sig1Orientation, StressGredient_Sig1Orientation, MinimumNormalStressAllowed, LoadingFaultCount,LoadingFaultInvert)
    

        FaultRakeAngle = Input_Bulk[:,10]
        P1, P2, P3, UnitVector_Normal, UnitVector_StrikeSlip, UnitVector_DipSlip, UnitVector_Slip = 
            RotVerts_UnitVectors(Input_Bulk, FaultCount, FaultRakeAngle) 
     
        FaultCenter = (P1+P2+P3)/3
        TractionVector = transpose(StressRatioXYZ * UnitVector_Normal')
        DotProduct = sum(UnitVector_Normal .* TractionVector, dims=2)
        Traction_Normal = DotProduct .* UnitVector_Normal
        Traction_Shear = TractionVector - Traction_Normal
        Vector_Horizontal = zeros(FaultCount,3)
        for ElemIdx = 1:FaultCount
            Vector_Horizontal[ElemIdx,:] = cross(UnitVector_Normal[ElemIdx,:], [0, 0, 1])

        end
        RakeAngle = acosd.(sum(Traction_Shear .* Vector_Horizontal, dims=2) ./ (norm.(eachrow(Traction_Shear)) .* norm.(eachrow(Vector_Horizontal)))) .+ 180
        Friction = norm.(eachrow(Traction_Shear)) ./ norm.(eachrow(Traction_Normal)) 

        if LoadingFaultCount > 0  && LoadingFaultInvert == 1
            RakeAngle[end-LoadingFaultCount+1:end] .= RakeAngle[end-LoadingFaultCount+1:end] .+ 180
        end

        for i in eachindex(Traction_Shear[:,3]); if Traction_Shear[i,3] < 0; RakeAngle[i] = 360 - RakeAngle[i] ; end; end
        for i in eachindex(RakeAngle); if RakeAngle[i] > 360; RakeAngle[i] -= 180; end; end
        
        for ElemIdx = 1:FaultCount
            # recalculate unitvectors with updated rake angle
            UnitVector_StrikeSlip[ElemIdx,:] = cross(UnitVector_Normal[ElemIdx,:], [0, 0, 1]) / 
                                                norm(cross(UnitVector_Normal[ElemIdx,:], [0, 0, 1]) )
            UnitVector_DipSlip[ElemIdx,:] = cross(UnitVector_StrikeSlip[ElemIdx,:], UnitVector_Normal[ElemIdx,:]) /
                                                norm(cross(UnitVector_StrikeSlip[ElemIdx,:], UnitVector_Normal[ElemIdx,:]))
            UnitVector_Slip[ElemIdx,:] = UnitVector_StrikeSlip[ElemIdx,:] * cosd(RakeAngle[ElemIdx]) + 
                                        UnitVector_DipSlip[ElemIdx,:] * sind(RakeAngle[ElemIdx])
        end
        Input_Bulk[:, 10] =  RakeAngle
        Input_Bulk[:, 16] =  Friction

        Input_Bulk[:,17] = norm.(eachrow(Traction_Normal)) .* StressOnSurface_Sig1Orientation
        Input_Bulk[:,18] = norm.(eachrow(Traction_Normal)) .* StressGredient_Sig1Orientation 
        TooSmallIdx = Input_Bulk[:,17] .< MinimumNormalStressAllowed
        Input_Bulk[TooSmallIdx,17] .= MinimumNormalStressAllowed


    return Input_Bulk, UnitVector_Normal, UnitVector_Slip, UnitVector_DipSlip, UnitVector_StrikeSlip
end


function CalculateFrictionStressHetero_T(Input_Bulk, MinFrictionAllowed,PrincipalStressRatioXYZ, 
                StressOnSurface_Sig1Orientation, StressGredient_Sig1Orientation, MinimumNormalStressAllowed,
                LoadingFaultCount, StressLocation, StressRotationStrike, StressRotationDip)
    

    
    FaultRakeAngle = Input_Bulk[:,10] # this will be over written after stresscalculation
    FaultCount = length(Input_Bulk[:,1])
    P1, P2, P3, UnitVector_Normal, UnitVector_StrikeSlip, UnitVector_DipSlip, UnitVector_Slip = 
        RotVerts_UnitVectors(Input_Bulk, FaultCount, FaultRakeAngle) ##### Here, Rake angle and Slip vector will be recalculated
    
    FaultCenter = (P1+P2+P3)/3
    for ElemIdx = 1:FaultCount

        UnitVector_Normal_i = UnitVector_Normal[ElemIdx,:]
        ReceiverPoint = FaultCenter[ElemIdx,:]
        
        ReceiverPrincipalStressRatioXYZ = zeros(3)
        ReceiverPrincipalStressRatioXYZ[1] = InverseDistWeighting(StressLocation, PrincipalStressRatioXYZ[:,1], ReceiverPoint, 1.0)
        ReceiverPrincipalStressRatioXYZ[2] = InverseDistWeighting(StressLocation, PrincipalStressRatioXYZ[:,2], ReceiverPoint, 1.0)
        ReceiverPrincipalStressRatioXYZ[3] = InverseDistWeighting(StressLocation, PrincipalStressRatioXYZ[:,3], ReceiverPoint, 1.0)
        ReceiverStressRotationStrike =  InverseDistWeighting(StressLocation, StressRotationStrike, ReceiverPoint, 1.0)
        ReceiverStressRotationDip =  InverseDistWeighting(StressLocation, StressRotationDip, ReceiverPoint, 1.0)
        PrincipalStressRatio = [-ReceiverPrincipalStressRatioXYZ[1] 0 0 
                                0 -ReceiverPrincipalStressRatioXYZ[2] 0
                                0 0 -ReceiverPrincipalStressRatioXYZ[3]]

        RotationMat_Strike =
        [cosd(ReceiverStressRotationStrike) -sind(ReceiverStressRotationStrike)  0
        sind(ReceiverStressRotationStrike) cosd(ReceiverStressRotationStrike) 0
        0  0  1]

        RotationMat_Dip =
        [1 0 0
        0 cosd(ReceiverStressRotationDip) -sind(ReceiverStressRotationDip) 
        0 sind(ReceiverStressRotationDip) cosd(ReceiverStressRotationDip)]

        RotationMat =  RotationMat_Strike * RotationMat_Dip
        StressRatioXYZ = RotationMat * PrincipalStressRatio * RotationMat'  

        
        TractionVector_i = StressRatioXYZ * UnitVector_Normal_i
        DotProduct_i = dot(UnitVector_Normal_i, TractionVector_i)            
        Traction_Normal_i = DotProduct_i .* UnitVector_Normal_i
        Traction_Shear_i = TractionVector_i - Traction_Normal_i
        Vector_Horizontal_i = cross(UnitVector_Normal_i, [0, 0, 1])
        RakeAngle = acosd(dot(Traction_Shear_i , Vector_Horizontal_i) / 
                    (norm(Traction_Shear_i ) * norm(Vector_Horizontal_i ))) + 180
        Friction = norm(Traction_Shear_i ) / norm(Traction_Normal_i ) 

        if Traction_Shear_i[3] < 0; RakeAngle = 360 - RakeAngle ; end
        if RakeAngle > 360; RakeAngle -= 180; end
        
        UnitVector_StrikeSlip[ElemIdx,:] = cross(UnitVector_Normal[ElemIdx,:], [0, 0, 1]) / 
                                            norm(cross(UnitVector_Normal[ElemIdx,:], [0, 0, 1]) )
        UnitVector_DipSlip[ElemIdx,:] = cross(UnitVector_StrikeSlip[ElemIdx,:], UnitVector_Normal[ElemIdx,:]) /
                                            norm(cross(UnitVector_StrikeSlip[ElemIdx,:], UnitVector_Normal[ElemIdx,:]))
        UnitVector_Slip[ElemIdx,:] = UnitVector_StrikeSlip[ElemIdx,:] * cosd(RakeAngle) + 
                                    UnitVector_DipSlip[ElemIdx,:] * sind(RakeAngle)


        Input_Bulk[ElemIdx, 10] =  RakeAngle
        Input_Bulk[ElemIdx, 16] =  Friction

        Input_Bulk[ElemIdx,17] = norm(Traction_Normal_i) * StressOnSurface_Sig1Orientation
        Input_Bulk[ElemIdx,18] = norm(Traction_Normal_i) * StressGredient_Sig1Orientation 
        if Input_Bulk[ElemIdx,17] < MinimumNormalStressAllowed; Input_Bulk[ElemIdx,17] = MinimumNormalStressAllowed; end

    end

    return Input_Bulk, UnitVector_Normal, UnitVector_Slip, UnitVector_DipSlip, UnitVector_StrikeSlip
end


function FindMu0_AdjV_R(Input_Bulk, MaximumTargetVelocity, V_p, V_r, ConstantMu0)
        MaxFric, MaxF_idx = findmax(Input_Bulk[:,14])
        V0 = 1e-9
        
        Fault_a = Input_Bulk[:,9] 
        Fault_b = Input_Bulk[:,10]
        Fault_Dc = Input_Bulk[:,11]
        Fault_Theta_i = Input_Bulk[:,12]
        if MaximumTargetVelocity > 0
            Mu0 = MaxFric - Fault_a[MaxF_idx] * log(MaximumTargetVelocity/V0) - Fault_b[MaxF_idx] * log(Fault_Theta_i[MaxF_idx]*V0 / Fault_Dc[MaxF_idx])
            println("Mu_0 Adjusted by target velocity. Mu0 = ", Mu0)
        else 
            Mu0 = ConstantMu0
            println("Input Mu_0 used. Mu0 =", Mu0)
        end
        Friction_0 = ones(length(Input_Bulk[:,9])) * Mu0
        Input_Bulk[:,13]  = V0 .* exp.( (Input_Bulk[:, 14] .- Friction_0 .- Fault_b .* log.(Fault_Theta_i .* V0./Fault_Dc)) ./ Fault_a)
        println("Velocity Adjusted. Mu_0 = ", Mu0, ". Maximum Velocity is ", maximum(Input_Bulk[:,13] ))
        if maximum(Input_Bulk[:,13]) > 1e-1
            println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            println("!!!!!!!!!!!! Warning: Maximum Velocity is larger than 0.1 m/s. Check the input parameters.!!!!!!!!!!")
            println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

        end
        PeakFriction = Mu0 + Fault_a[MaxF_idx] * log(V_p / 1e-9 ) + Fault_b[MaxF_idx] *log(1e-9 * Fault_Theta_i[MaxF_idx] ./ Fault_Dc[MaxF_idx])
        ResidualFriction = Mu0 + (Fault_a[MaxF_idx] - Fault_b[MaxF_idx]) * log.(V_r / 1e-9 )
    return Input_Bulk, Mu0, PeakFriction, ResidualFriction
end


function FindMu0_AdjV_T(Input_Bulk, MaximumTargetVelocity, V_p, V_r, ConstantMu0)
        MaxFric, MaxF_idx = findmax(Input_Bulk[:,16])
        V0 = 1e-9
        
        Fault_a = Input_Bulk[:,11] 
        Fault_b = Input_Bulk[:,12]
        Fault_Dc = Input_Bulk[:,13]
        Fault_Theta_i = Input_Bulk[:,14]
        Fault_friction = Input_Bulk[:,16]
        if MaximumTargetVelocity > 0
            Mu0 = MaxFric - Fault_a[MaxF_idx] * log(MaximumTargetVelocity/V0) - Fault_b[MaxF_idx] * log(Fault_Theta_i[MaxF_idx]*V0 / Fault_Dc[MaxF_idx])
            println("Mu_0 Adjusted by target velocity. Mu0 = ", Mu0)
        else 
            Mu0 = ConstantMu0
            println("Input Mu_0 used. Mu0 =", Mu0)
        end
        Friction_0 = ones(length(Input_Bulk[:,1])) * Mu0
        Input_Bulk[:,15]  = V0 .* exp.( (Fault_friction .- Friction_0 .- Fault_b .* log.(Fault_Theta_i .* V0./Fault_Dc)) ./ Fault_a)
        println("Velocity Adjusted. Mu_0 = ", Mu0, ". Maximum Velocity is ", maximum(Input_Bulk[:,15] ))
        if maximum(Input_Bulk[:,15]) > 1e-1
            println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            println("!!!!!!!!!!!! Warning: Maximum Velocity is larger than 0.1 m/s. Check the input parameters.!!!!!!!!!!")
            println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

        end
        PeakFriction = Mu0 + Fault_a[MaxF_idx] * log(V_p / 1e-9 ) + Fault_b[MaxF_idx] *log(1e-9 * Fault_Theta_i[MaxF_idx] ./ Fault_Dc[MaxF_idx])
        ResidualFriction = Mu0 + (Fault_a[MaxF_idx] - Fault_b[MaxF_idx]) * log.(V_r / 1e-9 )
    return Input_Bulk, Mu0, PeakFriction, ResidualFriction
end



function InverseDistWeighting(SourcePoint, SourceValue, ReceiverPoint, pValue)
    SourceCount = length(SourcePoint[:,1])
    Numerator = 0.0
    Denominator = 0.0
    for i=1:SourceCount
        DistanceInv = 1/ (norm(SourcePoint[i,:] - ReceiverPoint))^pValue
        if DistanceInv > 1e20
            DistanceInv = 1e20
        end
        Numerator = Numerator+ DistanceInv * SourceValue[i]
        Denominator = Denominator + DistanceInv
    end
    Value = Numerator/ Denominator
    return Value

end


function PlotPrincipalStresses(StressLocation, PrincipalStressRatioXYZ, StressRotationStrike, StressRotationDip, ax)

    StressCount = length(StressLocation[:,1])
    Linewidth = 2
    Arrow_length_ratio = 0.01
    PrinpalStressLength = 1e5
    for StressIdx=1:StressCount        
        UnlotatedVectorX =[PrincipalStressRatioXYZ[StressIdx,1], 0.0, 0.0]
        UnlotatedVectorY =[0.0, PrincipalStressRatioXYZ[StressIdx,2], 0.0]
        UnlotatedVectorZ =[0.0, 0.0, PrincipalStressRatioXYZ[StressIdx,3]]

        RotationMat_Strike=
        [cosd(StressRotationStrike[StressIdx]) -sind(StressRotationStrike[StressIdx])  0
        sind(StressRotationStrike[StressIdx]) cosd(StressRotationStrike[StressIdx]) 0
        0  0  1];

        RotationMat_Dip=
        [1 0 0
        0 cosd(StressRotationDip[StressIdx]) -sind(StressRotationDip[StressIdx])
        0 sind(StressRotationDip[StressIdx]) cosd(StressRotationDip[StressIdx])]


        RotatedStessX = RotationMat_Strike * RotationMat_Dip  * UnlotatedVectorX * PrinpalStressLength
        RotatedStessY = RotationMat_Strike * RotationMat_Dip  * UnlotatedVectorY * PrinpalStressLength
        RotatedStessZ = RotationMat_Strike * RotationMat_Dip  * UnlotatedVectorZ * PrinpalStressLength


        ax.quiver(StressLocation[StressIdx,1] + RotatedStessX[1], StressLocation[StressIdx,2]  + RotatedStessX[2], StressLocation[StressIdx,3] + RotatedStessX[3], 
        -RotatedStessX[1] , -RotatedStessX[2] , -RotatedStessX[3] ,
        color="k",arrow_length_ratio=Arrow_length_ratio, linewidth =Linewidth)

        ax.quiver(StressLocation[StressIdx,1] + RotatedStessY[1], StressLocation[StressIdx,2] + RotatedStessY[2],StressLocation[StressIdx,3]+ RotatedStessY[3], 
        -RotatedStessY[1] , -RotatedStessY[2] , -RotatedStessY[3] ,
        color="k",arrow_length_ratio=Arrow_length_ratio, linewidth =Linewidth)

        ax.quiver(StressLocation[StressIdx,1] + RotatedStessZ[1], StressLocation[StressIdx,2] + RotatedStessZ[2],StressLocation[StressIdx,3]+ RotatedStessZ[3], 
        -RotatedStessZ[1] , -RotatedStessZ[2] , -RotatedStessZ[3] ,
        color="k",arrow_length_ratio=Arrow_length_ratio, linewidth =Linewidth)

        ax.quiver(StressLocation[StressIdx,1] - RotatedStessX[1], StressLocation[StressIdx,2] - RotatedStessX[2], StressLocation[StressIdx,3] - RotatedStessX[3], 
        RotatedStessX[1] , RotatedStessX[2] , RotatedStessX[3] ,
        color="k",arrow_length_ratio=Arrow_length_ratio, linewidth =Linewidth)

        ax.quiver(StressLocation[StressIdx,1] - RotatedStessY[1], StressLocation[StressIdx,2] - RotatedStessY[2],StressLocation[StressIdx,3]- RotatedStessY[3], 
        RotatedStessY[1] , RotatedStessY[2] , RotatedStessY[3] ,
        color="k",arrow_length_ratio=Arrow_length_ratio, linewidth =Linewidth)

        ax.quiver(StressLocation[StressIdx,1] - RotatedStessZ[1], StressLocation[StressIdx,2] - RotatedStessZ[2],StressLocation[StressIdx,3] - RotatedStessZ[3], 
        RotatedStessZ[1] , RotatedStessZ[2] , RotatedStessZ[3] ,
        color="k",arrow_length_ratio=Arrow_length_ratio, linewidth =Linewidth)

        ax.set_aspect("equal")
    end
end


