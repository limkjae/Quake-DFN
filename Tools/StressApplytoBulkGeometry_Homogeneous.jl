
using DelimitedFiles
using PyPlot
using PyCall
using LinearAlgebra
using Statistics
pygui(true)
include("../Functions_BuildInputFile.jl")
include("../Results/Functions_Plot.jl")


InputBulkFileName="Input_BulkFaultGeometry.txt"

function ChangeBulk()

    ######################################## Inputs ########################################

    ###### build Principal Stress. Compression Positive. Only Ratio Matters! 
    PrincipalStressRatioX = 0.3
    PrincipalStressRatioY = 1.0
    PrincipalStressRatioZ = 0.6
    StressRotationStrike = -5 # degree
    StressRotationDip = 0 # degree

    MaximumTargetVelocity = 1e-11 # if this has value, the maximum velocity is set to this value. And Mu0 will be adjusted accordingly.
    ConstantTheta = 1e10 # if not zero, initial theta will be revised to this uniformly
    Fault_a_Rev = 0.0 # if not zero, RSF "a" value will be revised to this uniformly
    Fault_b_Rev = 0.0 # if not zero, RSF "b" value will be revised to this uniformly
    Fault_Dc_Rev = 0.0 # if not zero, RSF "Dc" value will be revised to this uniformly

    V_p = 1e-5 # When target velocity is set this will be used for peak friction plot
    V_r = 1e-2 # When target velocity is set this will be used for residual friction plot
    MinFrictionAllowed = 0.05
    
    MinimumNormalStressAllowed = 1e6
    StressOnSurface_Sig1Orientation = 10e6 # pascal
    StressGredient_Sig1Orientation = 0 # pascal/m


    FaultSegmentLength = 0 # if 0, segment length will be unchanged    


    LoadingFaultAdjust = 0 # if 0, Loading fault sense of slip will not be changed
    LoadingFaultInvert = 1 # if 1, loading fault sense of slip become inverted


    #####  FigureConfiguration  
    PlotPrincipalStress = 1 # 1:plot the principal stress (Aspect ratio: equal), 0: no
    PlotTractionVector = 0 # 1: plot traction, normal, shear vector, 0: no
    PlotLoadingFault = 0 # 1: plot loading fault, normal, shear vector, 0: no

    PlotRotation=[56,-109]
    Transparent = 1 # 1 for transparent fault plot
    Edge = 1 # 0 for no element boudary 
    MinMax_Axis=0

    StressVectorLocation = 0 # Autometically Adjusted when 0 
    PrinpalStressLength = 0 # Autometically Adjusted when 0 

    ########################################################################################

    MaxStressRatio =  maximum([PrincipalStressRatioX,PrincipalStressRatioY, PrincipalStressRatioZ])
    StressRotationStrike = StressRotationStrike + 0.001
    StressRotationDip = StressRotationDip + 0.001
    PrincipalStressRatioX = PrincipalStressRatioX / MaxStressRatio
    PrincipalStressRatioY = PrincipalStressRatioY / MaxStressRatio
    PrincipalStressRatioZ = PrincipalStressRatioZ / MaxStressRatio
    Sig1 = maximum([PrincipalStressRatioX, PrincipalStressRatioY, PrincipalStressRatioZ])
    Sig3 = minimum([PrincipalStressRatioX, PrincipalStressRatioY, PrincipalStressRatioZ])
    Input_Bulk=readdlm(InputBulkFileName)
    Switch_StrikeSlip_or_ReverseNormal = Input_Bulk[2,1] 
    Input_BulktoAdjust=Input_Bulk[4:end,:]
    BulkFaultCount = length(Input_BulktoAdjust[:,1])
    # Adjust LRRN to rake angle
    if Switch_StrikeSlip_or_ReverseNormal == 1
        for BulkIndex = 1: BulkFaultCount
            if Input_BulktoAdjust[BulkIndex,8] == -1.0
                Input_BulktoAdjust[BulkIndex,8] = 0.0
            else        
                Input_BulktoAdjust[BulkIndex,8] = 180.0
            end    
        end
    elseif Switch_StrikeSlip_or_ReverseNormal ==2
        for BulkIndex = 1: BulkFaultCount
            if Input_BulktoAdjust[BulkIndex,7] < 90.0
                if Input_BulktoAdjust[BulkIndex,8] == -1.0
                    Input_BulktoAdjust[BulkIndex,8] = 90.0
                else        
                    Input_BulktoAdjust[BulkIndex,8] = 270.0
                end    
            else 
                if Input_BulktoAdjust[BulkIndex,8] == -1.0
                    Input_BulktoAdjust[BulkIndex,8] = 270.0
                else        
                    Input_BulktoAdjust[BulkIndex,8] = 90.0
                end    
            end
        end
    end
        
    Input_BulktoAdjust = Input_BulktoAdjust[sortperm(Input_BulktoAdjust[:, 17], rev=false), :] # move the loading faults to the end
    if PlotLoadingFault == 0 
        LoadingFaultCountPlot = sum(x->x>0, Input_BulktoAdjust[:,17]) 
    else 
        LoadingFaultCountPlot = 0
    end

    if StressVectorLocation == 0
        StressVectorLocation = [mean(Input_BulktoAdjust[1:end - LoadingFaultCountPlot,1]), mean(Input_BulktoAdjust[1:end - LoadingFaultCountPlot,2]), 100]
    end 
    if PrinpalStressLength == 0
        PrinpalStressLength = maximum([maximum(Input_BulktoAdjust[1:end-LoadingFaultCountPlot,1]) - minimum(Input_BulktoAdjust[1:end-LoadingFaultCountPlot,1]), 
                                      maximum(Input_BulktoAdjust[1:end-LoadingFaultCountPlot,2]) - minimum(Input_BulktoAdjust[1:end-LoadingFaultCountPlot,2])])/2
    end

    #################### Calculate Stress for XYZ Frame #####################
    # Stres Negative for Compression
    PrincipalStressRatio = [-PrincipalStressRatioX 0 0 
                    0 -PrincipalStressRatioY 0
                    0 0 -PrincipalStressRatioZ]

    # Rotate the Stress along Maximum Stress Angle (XYZ corrdinate)
    RotationMat_Strike =
    [cosd(StressRotationStrike) -sind(StressRotationStrike)  0
    sind(StressRotationStrike) cosd(StressRotationStrike) 0
    0  0  1]

    RotationMat_Dip =
    [1 0 0
    0 cosd(StressRotationDip) -sind(StressRotationDip) 
    0 sind(StressRotationDip) cosd(StressRotationDip)]

    RotationMat =  RotationMat_Strike * RotationMat_Dip
    StressRatioXYZ = RotationMat * PrincipalStressRatio * RotationMat'  
    ########^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^########




    ############## Calcuate Stress and Frictions in Each Fault ###############
    RakeAngle = zeros(BulkFaultCount)
    ShearFric =  zeros(BulkFaultCount)
    NormalFric =  zeros(BulkFaultCount)
    NormStressSurface =  zeros(BulkFaultCount)
    NormStressGradient =  zeros(BulkFaultCount)
    for BulkIndex = 1:BulkFaultCount

        FaultStrikeAngle = Input_BulktoAdjust[BulkIndex, 6]
        FaultDipAngle = Input_BulktoAdjust[BulkIndex, 7]
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
        
        Input_BulktoAdjust[BulkIndex,15] = NormalFric[BulkIndex] .* StressOnSurface_Sig1Orientation
        if Input_BulktoAdjust[BulkIndex,15] <  MinimumNormalStressAllowed
            Input_BulktoAdjust[BulkIndex,15] = MinimumNormalStressAllowed
        end
        Input_BulktoAdjust[BulkIndex,16] = NormalFric[BulkIndex] .* StressGredient_Sig1Orientation 
        Friction = sqrt(StressOnFault[2,3]^2 + StressOnFault[1,3]^2) / abs(StressOnFault[3,3])

        if Friction < MinFrictionAllowed
            println("Fault ", BulkIndex, " has very small friction. Adjusted to Minimum Allowed Friction")
            ShearFric[BulkIndex] = ShearFric[BulkIndex]  * MinFrictionAllowed / Friction
            Friction = MinFrictionAllowed
        end
        if Input_BulktoAdjust[BulkIndex,17] > 0 && LoadingFaultInvert == 1
            Rake = Rake + 180
        end
        if Input_BulktoAdjust[BulkIndex,17] > 0 && LoadingFaultAdjust == 0
            # println("donotherin")
        else
            if StressOnFault[1,3] < 0;  Rake = Rake + 180;  end
            if Rake > 360; Rake = Rake - 360; end
            if Rake < 0; Rake = Rake + 360; end
            Input_BulktoAdjust[BulkIndex, 8] =  Rake
            Input_BulktoAdjust[BulkIndex, 14] =  Friction
        end
        # println(Rake, "   ", Friction)
    end

    if Fault_a_Rev != 0.0; Input_BulktoAdjust[:,9] .= Fault_a_Rev; end
    if Fault_b_Rev != 0.0; Input_BulktoAdjust[:,10] .= Fault_b_Rev; end
    if Fault_Dc_Rev != 0.0; Input_BulktoAdjust[:,11] .= Fault_Dc_Rev; end
    if ConstantTheta != 0.0; Input_BulktoAdjust[:,12] .= ConstantTheta; end

    PeakFriction = 0.0
    ResidualFriction = 0.0
    if MaximumTargetVelocity > 0 
        MaxFric, MaxF_idx = findmax(Input_BulktoAdjust[:,14])
        V0 = 1e-9
        
        Fault_a = Input_BulktoAdjust[:,9] 
        Fault_b = Input_BulktoAdjust[:,10]
        Fault_Dc = Input_BulktoAdjust[:,11]
        Fault_Theta_i = Input_BulktoAdjust[:,12]
        Mu0 = MaxFric - Fault_a[MaxF_idx] * log(MaximumTargetVelocity/V0) - Fault_b[MaxF_idx] * log(Fault_Theta_i[MaxF_idx]*V0 / Fault_Dc[MaxF_idx])
        Friction_0 = ones(length(Input_BulktoAdjust[:,9])) * Mu0
        Input_BulktoAdjust[:,13]  = V0 .* exp.( (Input_BulktoAdjust[:, 14] .- Friction_0 .- Fault_b .* log.(Fault_Theta_i .* V0./Fault_Dc)) ./ Fault_a)
        println("Velocity Adjusted. Mu_0 = ", Mu0, ". Maximum Velocity is ", maximum(Input_BulktoAdjust[:,13] ))
        PeakFriction = Mu0 + Fault_a[MaxF_idx] * log(V_p / 1e-9 ) + Fault_b[MaxF_idx] *log(1e-9 * Fault_Theta_i[MaxF_idx] ./ Fault_Dc[MaxF_idx])
        ResidualFriction = Mu0 + (Fault_a[MaxF_idx] - Fault_b[MaxF_idx]) * log.(V_r / 1e-9 )
    end

    if FaultSegmentLength > 0
        Input_BulktoAdjust[:,18]  .= FaultSegmentLength

    end
    SurvivedFaults = 0
    Input_BulktoAdjustFiltered = zeros(1,18) 
    for i=1:BulkFaultCount
        SurvivedFaults =+ 1
        if isnan(Input_BulktoAdjust[i,8])
            println("Slip sence of Bulk Element ",i, " cannot be defined with the given stress field. The element will be removed")
            println("This problem can be alleviated by adding anisotropy or small angles to the stress field applied")
        else
            Input_BulktoAdjustFiltered = [Input_BulktoAdjustFiltered; Input_BulktoAdjust[i,:]']
        end
    end
    Input_BulktoAdjustFiltered = Input_BulktoAdjustFiltered[2:end,:]
        ############################## Save File #############################

    open(InputBulkFileName, "w") do io
        write(io,"SwitchSS/RN\tShearMod\tPoissonRatio\tR_Density\tCrit_TooClose\tTooCloseNormal_Multiplier\tMinimum_NS\n")
        writedlm(io,[0.0   Input_Bulk[2,2]     Input_Bulk[2,3]      Input_Bulk[2,4]   Input_Bulk[2,5]      Input_Bulk[2,6]  MinimumNormalStressAllowed ])
        write(io, "Ctr_X\tCtr_Y\tCtr_Z\tSt_L\tDip_L\tStAng\tDipAng\tRake\ta\tb\tDc\tTheta_i\tV_i\tFric_i\tSig0\tSigGrad\tV_Const\tMaxLeng\n")
        writedlm(io, Input_BulktoAdjustFiltered)
    end

    figure(1)
    clf()
    PlotInput = Input_BulktoAdjust[1:end-LoadingFaultCountPlot,14]; ColorMinMax = 0
    MaxVaule, MinValue = FaultPlot_3D_Color_General(Input_BulktoAdjust[1:end-LoadingFaultCountPlot,1:3],
        Input_BulktoAdjust[1:end-LoadingFaultCountPlot,4], Input_BulktoAdjust[1:end-LoadingFaultCountPlot,5], Input_BulktoAdjust[1:end-LoadingFaultCountPlot,6], 
        Input_BulktoAdjust[1:end-LoadingFaultCountPlot,7], Input_BulktoAdjust[1:end-LoadingFaultCountPlot,8], PlotInput, 
        PlotRotation, MinMax_Axis, ColorMinMax, Transparent, Edge, 0)
        ax = subplot(projection="3d")
        xlabel("x")
        ylabel("y")
        plotforcbar=  scatter([1,1],[1,1],0.1, [MinValue,MaxVaule], cmap="jet")
        cbar  = colorbar(plotforcbar, pad=0.15)
        figure(1).canvas.draw()

    ax = PlotBulk_SenseOfSlip(0.0, Input_BulktoAdjust[1:end-LoadingFaultCountPlot,:], PlotRotation, Transparent, Edge, MinMax_Axis)
 
    MohrCircleCenter = Sig3 + (Sig1-Sig3)/2
    MohrCircleAngles = collect(0:5:180)
    MohrCircleX = MohrCircleCenter .+ cosd.(MohrCircleAngles) * (Sig1-Sig3)/2
    MohrCircleY = sind.(MohrCircleAngles)* (Sig1-Sig3)/2

    figure(3)
    clf()
    plot(MohrCircleX, MohrCircleY, color = [0.8, 0.8, 0.8])
    plot([0,1.2], [0,1.2] * 1.0, color = [0.8, 0.8, 0.8])
    plot([0,1.2], [0,1.2] * 0.8, color = [0.8, 0.8, 0.8])
    plot([0,1.2], [0,1.2] * 0.6, color = [0.8, 0.8, 0.8])
    plot([0,1.2], [0,1.2] * 0.4, color = [0.8, 0.8, 0.8])
    plot([0,1.2], [0,1.2] * 0.2, color = [0.8, 0.8, 0.8])
    plot([0,1.2], [0,1.2] * PeakFriction, color = "r")
    plot([0,1.2], [0,1.2] * ResidualFriction, color = "b")
    scatter(NormalFric, ShearFric,  facecolors="none", edgecolor="k")
    xlim([0,1.2])
    ylim([0,1.2])

    if PlotPrincipalStress ==1

        Linewidth = 2
        Arrow_length_ratio = 0.2

        UnlotatedVectorX =[PrincipalStressRatioX, 0.0, 0.0]
        UnlotatedVectorY =[0.0, PrincipalStressRatioY, 0.0]
        UnlotatedVectorZ =[0.0, 0.0, PrincipalStressRatioZ]


        RotationMat_Strike=
        [cosd(StressRotationStrike) -sind(StressRotationStrike)  0
        sind(StressRotationStrike) cosd(StressRotationStrike) 0
        0  0  1];

        RotationMat_Dip=
        [1 0 0
        0 cosd(StressRotationDip) -sind(StressRotationDip)
        0 sind(StressRotationDip) cosd(StressRotationDip)]

        RotatedStessX = RotationMat_Strike * RotationMat_Dip  * UnlotatedVectorX * PrinpalStressLength
        RotatedStessY = RotationMat_Strike * RotationMat_Dip  * UnlotatedVectorY * PrinpalStressLength
        RotatedStessZ = RotationMat_Strike * RotationMat_Dip  * UnlotatedVectorZ * PrinpalStressLength


        ax.quiver(StressVectorLocation[1] + RotatedStessX[1], StressVectorLocation[2] + RotatedStessX[2], StressVectorLocation[3]+ RotatedStessX[3], 
        -RotatedStessX[1] , -RotatedStessX[2] , -RotatedStessX[3] ,
        color="k",arrow_length_ratio=Arrow_length_ratio, linewidth =Linewidth)

        ax.quiver(StressVectorLocation[1] + RotatedStessY[1], StressVectorLocation[2] + RotatedStessY[2], StressVectorLocation[3]+ RotatedStessY[3], 
        -RotatedStessY[1] , -RotatedStessY[2] , -RotatedStessY[3] ,
        color="k",arrow_length_ratio=Arrow_length_ratio, linewidth =Linewidth)

        ax.quiver(StressVectorLocation[1] + RotatedStessZ[1], StressVectorLocation[2] + RotatedStessZ[2], StressVectorLocation[3]+ RotatedStessZ[3], 
        -RotatedStessZ[1] , -RotatedStessZ[2] , -RotatedStessZ[3] ,
        color="k",arrow_length_ratio=Arrow_length_ratio, linewidth =Linewidth)

        ax.quiver(StressVectorLocation[1] - RotatedStessX[1], StressVectorLocation[2] - RotatedStessX[2], StressVectorLocation[3] - RotatedStessX[3], 
        RotatedStessX[1] , RotatedStessX[2] , RotatedStessX[3] ,
        color="k",arrow_length_ratio=Arrow_length_ratio, linewidth =Linewidth)

        ax.quiver(StressVectorLocation[1] - RotatedStessY[1], StressVectorLocation[2] - RotatedStessY[2], StressVectorLocation[3] - RotatedStessY[3], 
        RotatedStessY[1] , RotatedStessY[2] , RotatedStessY[3] ,
        color="k",arrow_length_ratio=Arrow_length_ratio, linewidth =Linewidth)

        ax.quiver(StressVectorLocation[1] - RotatedStessZ[1], StressVectorLocation[2] - RotatedStessZ[2], StressVectorLocation[3] - RotatedStessZ[3], 
        RotatedStessZ[1] , RotatedStessZ[2] , RotatedStessZ[3] ,
        color="k",arrow_length_ratio=Arrow_length_ratio, linewidth =Linewidth)
        
        ax.set_aspect("auto")
    end

    ################ Traction Vector Plots ###############
    if PlotTractionVector == 1
        LineLength = vec(minimum([Input_BulktoAdjust[:,4] Input_BulktoAdjust[:,5]], dims=2)./2)
        for BulkIndex =1:BulkFaultCount     

            FaultStrikeAngle = Input_BulktoAdjust[BulkIndex, 6]
            FaultDipAngle = Input_BulktoAdjust[BulkIndex, 7]
            ############# Calcualte Stress from Traction Vector ###########
                # println(FaultStrikeAngle[BulkIndex])
            VectorRotation_Strike = 
            [cosd(FaultStrikeAngle) -sind(FaultStrikeAngle)  0
            sind(FaultStrikeAngle) cosd(FaultStrikeAngle) 0
            0  0  1]
            VectorRotation_Dip = 
            [1 0 0
            0 cosd(FaultDipAngle) -sind(FaultDipAngle) 
            0 sind(FaultDipAngle) cosd(FaultDipAngle)]
            NormalVector =  VectorRotation_Strike * VectorRotation_Dip* [0, 0, -1]
            TractionVector = StressRatioXYZ * NormalVector
            NormalMagnitude = dot(TractionVector, NormalVector)
            ShearVector = TractionVector - NormalVector * NormalMagnitude
            ShearMagnitude = sqrt(norm(TractionVector)^2 -NormalMagnitude^2)
            # Friction = ShearMagnitude/NormalMagnitude

            ax.quiver(Input_BulktoAdjust[BulkIndex,1], Input_BulktoAdjust[BulkIndex,2], -Input_BulktoAdjust[BulkIndex,3], 
            TractionVector[1] * LineLength[BulkIndex], TractionVector[2] * LineLength[BulkIndex], TractionVector[3] * LineLength[BulkIndex],
                color="r",arrow_length_ratio=0.2)
                
            ax.quiver(Input_BulktoAdjust[BulkIndex,1], Input_BulktoAdjust[BulkIndex,2], -Input_BulktoAdjust[BulkIndex,3], 
            NormalVector[1] * LineLength[BulkIndex], NormalVector[2] * LineLength[BulkIndex], NormalVector[3] * LineLength[BulkIndex],
                color="g",arrow_length_ratio=0.2)
                
            ax.quiver(Input_BulktoAdjust[BulkIndex,1], Input_BulktoAdjust[BulkIndex,2], -Input_BulktoAdjust[BulkIndex,3], 
            ShearVector[1] * LineLength[BulkIndex], ShearVector[2] * LineLength[BulkIndex], ShearVector[3] * LineLength[BulkIndex],
                color="b",arrow_length_ratio=0.2)

        end
    end


end

ChangeBulk()