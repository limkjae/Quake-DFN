
using DelimitedFiles
using PyPlot
using PyCall
pygui(true)
include("../Functions_BuildInputFile.jl")
include("../Results/Functions_Plot.jl")


InputBulkFileName="Input_BulkFaultGeometry.txt"


# build Principal Stress. Compression Positive. Only Ratio Matters!
PrincipalStressRatioX = 0.2
PrincipalStressRatioY = 0.2
PrincipalStressRatioZ = 1.0
StressRotationStrike = 0 # degree
StressRotationDip = 45 # degree

PlotPrincipalStress = 1 # 1:plot the principal stress (Aspect ratio: equal), 0: no
PrinpalStressLength = 500 # if plotted, maximum length is...
StressVectorLocation = [0,0,0]

# MaxStressRatio = maximum([PrincipalStressRatioX , PrincipalStressRatioY, PrincipalStressRatioZ])
# PrincipalStressRatioX , PrincipalStressRatioY, PrincipalStressRatioZ = 
# PrincipalStressRatioX/MaxStressRatio , PrincipalStressRatioY/MaxStressRatio, PrincipalStressRatioZ /MaxStressRatio

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


Input_Bulk=readdlm(InputBulkFileName)
Input_BulktoAdjust=Input_Bulk[4:end,:]
BulkFaultCount = length(Input_BulktoAdjust[:,1])

        ##  1.Ctr_X     2.Ctr_Y 3.Ctr_Z 4.St_L	    5.Dip_L	    6.StAng	    7.DipAng	8.Rake
        ##  9.a         10.b	11.Dc	12.Theta_i	13. V_i     14. Friction_i 15.NormalStress at surface [Pa]  
        ##  16. NoarmalStress Gradient [Pa] 17. V_Const     18. Minimum Segment Length
RakeAngle = zeros(BulkFaultCount)

for BulkIndex = 1:BulkFaultCount
    # BulkIndex = 17
    # BulkIndex = 1
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
    println(StressOnFault)
    if StressOnFault[1,3] < 0;  Rake = Rake + 180;  end
    if Rake > 360; Rake = Rake - 360; end
    if Rake < 0; Rake = Rake + 360; end
    Input_BulktoAdjust[BulkIndex, 8] =  Rake
    Input_BulktoAdjust[BulkIndex, 14] =  Friction
end

# println(StressOnFault)
# println(Rake)


    ############################## Save File #############################

open(InputBulkFileName, "w") do io
    write(io,"SwitchSS/RN\tShearMod\tPoissonRatio\tR_Density\tCrit_TooClose\tTooCloseNormal_Multiplier\tMinimum_NS\n")
    writedlm(io,[0.0   Input_Bulk[2,2]     Input_Bulk[2,3]      Input_Bulk[2,4]   Input_Bulk[2,5]      Input_Bulk[2,6]  Input_Bulk[2,7] ])
    write(io, "Ctr_X\tCtr_Y\tCtr_Z\tSt_L\tDip_L\tStAng\tDipAng\tRakeAngle\ta\tb\tDc\tTheta_i\tV_i\tFric_i\tSig0\tSigGrad\tV_Const\tMaxLeng\n")
    writedlm(io, Input_BulktoAdjust)
end




include("../Plot_BulkFaultGeometry.jl")



if PlotPrincipalStress ==1

    Linewidth = 3
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