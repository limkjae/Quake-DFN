
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

    ###### build Principal Stress. Compression Positive. Only Ratio Matters!
    PrincipalStressRatioX = 0.3
    PrincipalStressRatioY = 1.0
    PrincipalStressRatioZ = 0.5
    StressRotationStrike = -10 # degree
    StressRotationDip = 5 # degree

    ##############  FigureConfiguration  ##############

    PlotPrincipalStress = 1 # 1:plot the principal stress (Aspect ratio: equal), 0: no
    PlotTractionVector = 0 # 1: plot traction, normal, shear vector, 0: no

    PlotRotation=[30,-50]
    Transparent = 1 # 1 for transparent fault plot
    Edge = 1 # 0 for no element boudary 
    MinMax_Axis=0

    StressVectorLocation = 0 # Autometically Adjusted when 0 
    PrinpalStressLength = 0 # Autometically Adjusted when 0 
    ###################################################


    Input_Bulk=readdlm(InputBulkFileName)
    Input_BulktoAdjust=Input_Bulk[4:end,:]
    Input_BulktoAdjust = Input_BulktoAdjust[sortperm(Input_BulktoAdjust[:, 16], rev=true), :] # move the loading faults to the end
    BulkFaultCount = length(Input_BulktoAdjust[:,1])


    if StressVectorLocation == 0
        StressVectorLocation = [mean(Input_BulktoAdjust[:,1]), mean(Input_BulktoAdjust[:,2]), 100]
    end 
    if PrinpalStressLength == 0
        PrinpalStressLength = maximum([maximum(Input_BulktoAdjust[:,1]) - minimum(Input_BulktoAdjust[:,1]), 
                                      maximum(Input_BulktoAdjust[:,2]) - minimum(Input_BulktoAdjust[:,2])])/5

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
        Friction = sqrt(StressOnFault[2,3]^2 + StressOnFault[1,3]^2) / abs(StressOnFault[3,3])
        println(Rake, "   ", Friction)
        # println(StressOnFault)
        if StressOnFault[1,3] < 0;  Rake = Rake + 180;  end
        if Rake > 360; Rake = Rake - 360; end
        if Rake < 0; Rake = Rake + 360; end
        Input_BulktoAdjust[BulkIndex, 8] =  Rake
        Input_BulktoAdjust[BulkIndex, 14] =  Friction

    end

        ############################## Save File #############################

    open(InputBulkFileName, "w") do io
        write(io,"SwitchSS/RN\tShearMod\tPoissonRatio\tR_Density\tCrit_TooClose\tTooCloseNormal_Multiplier\tMinimum_NS\n")
        writedlm(io,[0.0   Input_Bulk[2,2]     Input_Bulk[2,3]      Input_Bulk[2,4]   Input_Bulk[2,5]      Input_Bulk[2,6]  Input_Bulk[2,7] ])
        write(io, "Ctr_X\tCtr_Y\tCtr_Z\tSt_L\tDip_L\tStAng\tDipAng\tRakeAngle\ta\tb\tDc\tTheta_i\tV_i\tFric_i\tSig0\tSigGrad\tV_Const\tMaxLeng\n")
        writedlm(io, Input_BulktoAdjust)
    end


    figure(1)
    clf()
    PlotInput = Input_BulktoAdjust[:,14]; ColorMinMax = 0
    MaxVaule, MinValue = FaultPlot_3D_Color_General(Input_BulktoAdjust[:,1:3],
        Input_BulktoAdjust[:,4], Input_BulktoAdjust[:,5], Input_BulktoAdjust[:,6], Input_BulktoAdjust[:,7], Input_BulktoAdjust[:,8], PlotInput, 
        PlotRotation, MinMax_Axis, ColorMinMax, Transparent, Edge, 0)
        ax = subplot(projection="3d")
        xlabel("x")
        ylabel("y")
        plotforcbar=  scatter([1,1],[1,1],0.1, [MinValue,MaxVaule], cmap="jet")
        cbar  = colorbar(plotforcbar, pad=0.15)
        figure(1).canvas.draw()

    ax = PlotBulk_SenseOfSlip(0.0, Input_BulktoAdjust, PlotRotation, Transparent, Edge, MinMax_Axis)
 

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
        
        ax.set_aspect("equal")
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