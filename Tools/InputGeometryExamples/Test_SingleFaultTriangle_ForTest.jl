
using LinearAlgebra
using Statistics
using PyPlot
using PyCall
using JLD2
using PyPlot
using PyCall
using DelimitedFiles


pygui(true)
# include("Functions_TDstressHS.jl")

OutputFileName="Input_BulkFaultGeometry.txt"

TotalElemCount = 2
FaultLength = 100
Rake = 0
a = 0.003
b = 0.006
Dc = 3e-4
Thata_i = 1e5
V_i = 1e-15
Friction_i = 0.6
NormalStressAtSurface = 10e6
NormalStressGradient = 0
V_Const = 0.0
MinimumSegmentLength = 100000.0


Rake_all = ones(TotalElemCount) * Rake
a_all = ones(TotalElemCount) * a 
b_all = ones(TotalElemCount) * b
Dc_all = ones(TotalElemCount) * Dc
Thata_i_all = ones(TotalElemCount) * Thata_i
V_i_all = ones(TotalElemCount) * V_i
Friction_i_all = ones(TotalElemCount) * Friction_i
NormalStressAtSurface_all = ones(TotalElemCount) * NormalStressAtSurface
NormalStressGradient_all = ones(TotalElemCount) * NormalStressGradient
V_Const_all = ones(TotalElemCount) * V_Const
    
# FourVerts = [-FaultLength 0 0; -FaultLength 0 -FaultLength; FaultLength 0.0 -FaultLength; FaultLength 0.0 0]
# FourVerts = [ 0 -FaultLength 0;  0 -FaultLength -FaultLength; 0.0  FaultLength -FaultLength;  0.0 FaultLength 0]
# FourVerts = [-FaultLength -FaultLength 0;  -FaultLength  -FaultLength -FaultLength; 0.0  FaultLength -FaultLength;  0.0 FaultLength 0]
# FourVerts = [FaultLength -FaultLength 0;  FaultLength  -FaultLength -FaultLength; 0.0  FaultLength -FaultLength;  0.0 FaultLength 0]
FourVerts = [FaultLength FaultLength*5 -FaultLength*5;  FaultLength  -FaultLength -FaultLength; 0.0  FaultLength -FaultLength;  0.0 FaultLength 0]
LoadingV = 1e-9
LPInput1 = [FourVerts[1,:]; FourVerts[2,:]; FourVerts[3,:]; Rake+180;  a; b; Dc; Thata_i;  V_i;  Friction_i;  NormalStressAtSurface;  NormalStressGradient;  0.0; MinimumSegmentLength]
LPInput2 = [FourVerts[1,:]; FourVerts[3,:]; FourVerts[4,:]; Rake+180;  a; b; Dc; Thata_i;  V_i;  Friction_i;  NormalStressAtSurface;  NormalStressGradient;  0.0; MinimumSegmentLength]
InputFile = [LPInput1'; LPInput2']

P1 = [InputFile[:,1] InputFile[:,2] InputFile[:,3]]
P2 = [InputFile[:,4] InputFile[:,5] InputFile[:,6]]
P3 = [InputFile[:,7] InputFile[:,8] InputFile[:,9]]

InputProperty = InputFile[:,3]
MinValue = minimum(InputProperty)
MaxValue = maximum(InputProperty)
figure(1)
fig = figure(1)
clf()
art3d = PyObject(PyPlot.art3D)
ax = subplot(projection="3d")
for ElemIdx = 1:TotalElemCount
    cm = get_cmap(:jet)
    PlotValue=(InputProperty[ElemIdx]-MinValue)/(MaxValue-MinValue)

    face_color = [cm(PlotValue)[1], cm(PlotValue)[2],cm(PlotValue)[3], 0.1]

    verts = ((P1[ElemIdx,:],P2[ElemIdx,:],P3[ElemIdx,:]), )
    p3c = PyObject(art3d.Poly3DCollection(verts))
    pycall(ax.add_collection3d, PyAny, p3c)

    # face_color = [0.3, 0.8, 0.3, 0.5]         
    edge_color = [0.2, 0.2, 0.2, 1.0]

    pycall(p3c.set_facecolor, PyAny, face_color)
    pycall(p3c.set_edgecolor, PyAny, edge_color)
    ax.view_init(45, -30)

end
xlabel("x")
ylabel("y")

ax.set_aspect("equal")

SwitchSSRN = 0
ShearMod = 20e9
PoissonRatio = 0.25
R_Density=2670.0
Crit_TooClose= 1.05
TooCloseNormal_Multiplier = 0.6
MinNormalStress = 2e6


    open(OutputFileName, "w") do io
        write(io,"SwitchSS/RN\tShearMod\tPoissonRatio\tR_Density\tCrit_TooClose\tTooCloseNormal_Multiplier\tMinimum_NS\n")
        writedlm(io,[SwitchSSRN   ShearMod    PoissonRatio     R_Density  Crit_TooClose     TooCloseNormal_Multiplier MinNormalStress])
        write(io, "P1X\tP1Y\tP1Z\tP2X\tP2Y\tP2Z\tP3X\tP3Y\tP3Z\tRakeAngle\ta\tb\tDc\tTheta_i\tV_i\tFric_i\tSig0\tSigGrad\tV_Const\tMaxLeng\n")
        writedlm(io, InputFile)
    end;

    println("Saved File Name: ",OutputFileName)



# CoM = (P1+P2+P3)/3
# InputProperty = CoM[:,3] 
# MaxValue=maximum(InputProperty)
# MinValue=minimum(InputProperty)




# # SpecificStressStiffness = zeros(TotalElemCount,TotalElemCount)
# println("Compiling Stiffness Matrix Function. This may take time if first run")

# for ElemIdx = 1:TotalElemCount

#     Stress,Strain = TDstressHS(FaultCenter[:,1],FaultCenter[:,2],FaultCenter[:,3],P1[ElemIdx,:],P2[ElemIdx,:],P3[ElemIdx,:],Ss,Ds,Ts,mu,lambda)
#     println(ElemIdx)   
#     for ElemIdx2 = 1:TotalElemCount
#         Stress_i = [Stress[ElemIdx2,1] Stress[ElemIdx2,4] Stress[ElemIdx2,5]
#                     Stress[ElemIdx2,4] Stress[ElemIdx2,2] Stress[ElemIdx2,6]
#                     Stress[ElemIdx2,5] Stress[ElemIdx2,6] Stress[ElemIdx2,3]]

#         TVector = Stress_i * UnitVector_Normal[ElemIdx2,:]
#         Stress_Normal = dot(TVector, UnitVector_Normal[ElemIdx2,:])
#         Stress_SS = dot(TVector, UnitVector_SS[ElemIdx2,:])
#         Stress_Dip = dot(TVector, UnitVector_Dip[ElemIdx2,:])
#         Stress_Shear = dot(TVector, UnitVector_Slip[ElemIdx2,:])
#         StiffnessMatrix_Normal[ElemIdx2, ElemIdx] = Stress_Normal
#         StiffnessMatrix_Shear[ElemIdx2, ElemIdx] = Stress_Shear
#         # SpecificStressStiffness[ElemIdx, ElemIdx2] =  Stress[ElemIdx2,6]
#     end
# end


# OutputFileName="Verification/StiffnessCurvedTriangleSSFiner.jld2"
# save(OutputFileName, 
# "StiffnessMatrixShear", StiffnessMatrix_Shear, "StiffnessMatrixNormal", StiffnessMatrix_Normal, 
# "P1",P1,"P2", P2, "P3",P3)
# println("Saved File Name: ",OutputFileName)
