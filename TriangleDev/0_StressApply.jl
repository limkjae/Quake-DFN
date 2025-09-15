
using LinearAlgebra
using Statistics
using PyPlot
using PyCall
using JLD2
using PyPlot
using PyCall
pygui(true)
using DelimitedFiles
include("../Functions_BuildInputFile.jl")


# include("Functions_TDstressHS.jl")
InputBulkFileName="Input_BulkFaultGeometry.txt"


###### build Principal Stress. Compression Positive. Only Ratio Matters! 
PrincipalStressRatioX = 0.3
PrincipalStressRatioY = 1.0
PrincipalStressRatioZ = 0.3
StressRotationStrike = 0 # degree
StressRotationDip = 60 # degree

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



##### Processing Principal Stresses
MaxStressRatio =  maximum([PrincipalStressRatioX,PrincipalStressRatioY, PrincipalStressRatioZ])
StressRotationStrike = StressRotationStrike + 0.001
StressRotationDip = StressRotationDip + 0.001
PrincipalStressRatioX = PrincipalStressRatioX / MaxStressRatio
PrincipalStressRatioY = PrincipalStressRatioY / MaxStressRatio
PrincipalStressRatioZ = PrincipalStressRatioZ / MaxStressRatio
Sig1 = maximum([PrincipalStressRatioX, PrincipalStressRatioY, PrincipalStressRatioZ])
Sig3 = minimum([PrincipalStressRatioX, PrincipalStressRatioY, PrincipalStressRatioZ])

#################### Calculate Stress for XYZ Frame #####################
# Stres Negative for Compression
PrincipalStressRatio = [-PrincipalStressRatioX 0 0 
                        0 -PrincipalStressRatioY 0
                        0 0 -PrincipalStressRatioZ]
StressRotationStrike = -StressRotationStrike
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
StressRatioCartesian = RotationMat * PrincipalStressRatio * RotationMat' 
########^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^########


Input_Bulk=readdlm(InputBulkFileName)
if size(Input_Bulk, 2) == 18
    println("Rectangle")
    RorT = "R"
elseif  size(Input_Bulk, 2) == 20
    println("Triangle")
    RorT = "T"
else 
    error("Input Bulk Fault Geometry file should have 18 or 20 columns")
end

########################## Read and Remove Header then Segmentize #########################
SwitchSSRN = Input_Bulk[2,1] 
ShearMod = Input_Bulk[2,2]
PoissonRatio = Input_Bulk[2,3]
R_Density = Input_Bulk[2,4]
Crit_TooClose= Input_Bulk[2,5]
TooCloseNormal_Multiplier= Input_Bulk[2,6]
MinNormalStress=Input_Bulk[2,7]



Input_Segment, LoadingFaultCount, ShearModulus, PoissonRatio, RockDensity, 
Switch_StrikeSlip_or_ReverseNormal, DropCrit, DropCritNormalStressMultiplier, MinimumNS, RorT, FaultCount =
    ReadBulkInput(InputBulkFileName)

FaultCenter, Fault_a, Fault_b, Fault_Dc, Fault_Theta_i, Fault_V_i, Fault_Friction_i, Fault_NormalStress, 
Fault_V_Const, Fault_BulkIndex, FaultLengthStrike, FaultLengthDip, FaultStrikeAngle, 
FaultDipAngle, FaultRakeAngle, FaultLengthStrike_Bulk, FaultLengthDip_Bulk, NormalStiffnessZero = 
    ReadSegmentInput(Input_Segment, FaultCount, RorT) 

P1, P2, P3, UnitVector_Normal, UnitVector_StrikeSlip, UnitVector_DipSlip, UnitVector_Slip = 
    RotVerts_UnitVectors(Input_Segment, FaultCount, FaultRakeAngle) 

FaultCenter = (P1+P2+P3)/3

TractionVector = transpose(StressRatioCartesian * UnitVector_Normal')
DotProduct = sum(UnitVector_Normal .* TractionVector, dims=2)
Traction_Normal = DotProduct .* UnitVector_Normal
Traction_Shear = TractionVector - Traction_Normal
UnitVector_Horizontal = zeros(FaultCount,3)
for ElemIdx = 1:FaultCount
    UnitVector_Horizontal[ElemIdx,:] = cross(UnitVector_Normal[ElemIdx,:], [0, 0, 1])
end

RakeAngle = acosd.(sum(Traction_Shear .* UnitVector_Horizontal, dims=2) ./ (norm.(eachrow(Traction_Shear)) .* norm.(eachrow(UnitVector_Horizontal)))) .+ 180

Friction = norm.(eachrow(Traction_Shear)) ./ norm.(eachrow(Traction_Normal)) 
if LoadingFaultCount > 0
    RakeAngle[end-LoadingFaultCount+1:end] .= RakeAngle[end-LoadingFaultCount+1:end] .+ 180
end
for i in eachindex(Traction_Shear[:,3]); if Traction_Shear[i,3] < 0; RakeAngle[i] = 360 - RakeAngle[i] ; end; end
for i in eachindex(RakeAngle); if RakeAngle[i] > 360; RakeAngle[i] -= 180; end; end

figure(10); clf(); plot(RakeAngle, ".")

InputProperty = Friction
# InputProperty = zeros(FaultCount)
MaxValue=maximum(InputProperty)
MinValue=minimum(InputProperty)

 

ArrowLength = 500
figure(2)
fig = figure(2)
clf()
art3d = PyObject(PyPlot.art3D)
ax = subplot(projection="3d")
for ElemIdx = 1:FaultCount - LoadingFaultCount
    cm = get_cmap(:jet)
    PlotValue=(InputProperty[ElemIdx]-MinValue)/(MaxValue-MinValue)

    face_color = [cm(PlotValue)[1], cm(PlotValue)[2],cm(PlotValue)[3], 0.5]

    verts = ((P1[ElemIdx,:],P2[ElemIdx,:],P3[ElemIdx,:]), )
    p3c = PyObject(art3d.Poly3DCollection(verts))
    pycall(ax.add_collection3d, PyAny, p3c)

    # face_color = [0.3, 0.8, 0.3, 0.5]         
    edge_color = [0.2, 0.2, 0.2, 1.0]

    pycall(p3c.set_facecolor, PyAny, face_color)
    pycall(p3c.set_edgecolor, PyAny, edge_color)
    ax.view_init(45, -30)

    # ### NormalVector Plot
    # ax.quiver(FaultCenter[ElemIdx,1], FaultCenter[ElemIdx,2], FaultCenter[ElemIdx,3], 
    #     UnitVector_Normal[ElemIdx,1] * ArrowLength, UnitVector_Normal[ElemIdx,2] * ArrowLength,
    #     UnitVector_Normal[ElemIdx,3] * ArrowLength,
    #     color="k",arrow_length_ratio=0.2)

        
    ### Horozontal Vector Plot
    ax.quiver(FaultCenter[ElemIdx,1], FaultCenter[ElemIdx,2], FaultCenter[ElemIdx,3], 
        Traction_Shear[ElemIdx,1] * ArrowLength, Traction_Shear[ElemIdx,2] * ArrowLength,
        Traction_Shear[ElemIdx,3] * ArrowLength,
        color="r",arrow_length_ratio=0.2)
end
xlabel("x")
ylabel("y")

StressVectorLocation = [mean(FaultCenter[:,1]), mean(FaultCenter[:,2]), 100]

PrinpalStressLength = 500
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


Input_Segment[:,10] = RakeAngle
Input_Segment[:,16] = Friction

open(InputBulkFileName, "w") do io
    write(io,"SwitchSS/RN\tShearMod\tPoissonRatio\tR_Density\tCrit_TooClose\tTooCloseNormal_Multiplier\tMinimum_NS\n")
    writedlm(io,[SwitchSSRN   ShearMod    PoissonRatio     R_Density  Crit_TooClose     TooCloseNormal_Multiplier MinNormalStress])
    write(io, "P1X\tP1Y\tP1Z\tP2X\tP2Y\tP2Z\tP3X\tP3Y\tP3Z\tRakeAngle\ta\tb\tDc\tTheta_i\tV_i\tFric_i\tSig0\tSigGrad\tV_Const\tMaxLeng\n")
    writedlm(io, Input_Segment)
end
