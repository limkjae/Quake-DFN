
using LinearAlgebra
using Statistics
using PyPlot
using PyCall
using JLD2
using PyPlot
using PyCall


pygui(true)
# include("Functions_TDstressHS.jl")


Length = 1000
Depth = 1000
Amplitude = 300
XElemCount = 50
ZElemCount = 50

Ss = 1
Ds = 0
Ts = 0

mu = 20e9
lambda = 20e9

    ############################# Write Bulk Input #################################
    ######++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++######
    #### Bulk File Order 
    ####  1.Ctr_X     2.Ctr_Y 3.Ctr_Z 4.St_L	    5.Dip_L	    6.StAng	    7.DipAng	8.LR/RN (-1: LL / +1 RL or -1:Reverse / +1 Normal)
    ####  9.a         10.b	11.Dc	12.Theta_i	13. V_i     14. Friction_i 15.NormalStress at surface [Pa]  
    ####  16. NoarmalStress Gradient [Pa] 17. V_Const     18. Minimum Segment Length
    ######++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++######
    



XElemLength = Length / XElemCount
ZElemLength = Depth / ZElemCount
TotalElemCount = XElemCount * ZElemCount * 2
FaultCenter = zeros(TotalElemCount,3)
UnitVector_Normal = zeros(TotalElemCount,3)
UnitVector_SS = zeros(TotalElemCount,3)
UnitVector_Dip = zeros(TotalElemCount,3)
UnitVector_Slip = zeros(TotalElemCount,3)

StiffnessMatrix_Shear = zeros(TotalElemCount,TotalElemCount)
StiffnessMatrix_Normal= zeros(TotalElemCount,TotalElemCount)

P1 = zeros(TotalElemCount,3)
P2 = zeros(TotalElemCount,3)
P3 = zeros(TotalElemCount,3)

for ElemIdx = 1:XElemCount * ZElemCount * 2

    ElemRow = ceil(ElemIdx / XElemCount /2)
    OrderInRow = ElemIdx - (ElemRow-1) * XElemCount* 2
    BlockInRow = floor((OrderInRow-1) / 2) + 1

    Xvert1 = (BlockInRow -1) * XElemLength
    Xvert2 = BlockInRow * XElemLength
    Zvert1 = -(ElemRow -1) * ZElemLength
    Zvert2 = - ElemRow * ZElemLength 
    
    if mod(OrderInRow,2) == 1
        x = [Xvert1 Xvert1 Xvert2]
        z = [Zvert1 Zvert2 Zvert2]

    else
        x = [Xvert1 Xvert2 Xvert2]
        z = [Zvert1 Zvert1 Zvert2]
    end    
    y = - sind.(x / Length * 180) * Amplitude
    

    FaultCenter[ElemIdx,:] = [mean(x), mean(y), mean(z)] 
    
    P1_i = [x[1], y[1], z[1]]
    P2_i = [x[2], y[2], z[2]]
    P3_i = [x[3], y[3], z[3]]
    UnitVector_Normal_i = cross(P2_i-P1_i, P3_i-P1_i)
    
    if angle(UnitVector_Normal_i[1] + UnitVector_Normal_i[2]*im) <= 0 
        P_temp = P1_i
        P1_i = P2_i
        P2_i = P_temp        
    end

    P1[ElemIdx,:] = P1_i
    P2[ElemIdx,:] = P2_i
    P3[ElemIdx,:] = P3_i
    
    UnitVector_Normal_i = cross(P2_i-P1_i, P3_i-P1_i)
    UnitVector_Normal_i = UnitVector_Normal_i/norm(UnitVector_Normal_i)
    UnitVector_SS_i = cross(UnitVector_Normal_i, [0,0,1])
    UnitVector_Dip_i = cross(UnitVector_Normal_i, UnitVector_SS_i)
    
    UnitVector_Normal[ElemIdx,:] = UnitVector_Normal_i
    UnitVector_SS[ElemIdx,:] = UnitVector_SS_i
    UnitVector_Dip[ElemIdx,:] = UnitVector_Dip_i

    UnitVector_Slip[ElemIdx,:] = (UnitVector_SS_i * Ss + UnitVector_Dip_i * Ds) / norm(UnitVector_SS_i * Ss + UnitVector_Dip_i * Ds)

end
CoM = (P1+P2+P3)/3
InputProperty = CoM[:,3] 
MaxValue=maximum(InputProperty)
MinValue=minimum(InputProperty)

figure(1)
fig = figure(1)
clf()
art3d = PyObject(PyPlot.art3D)
ax = subplot(projection="3d")
for ElemIdx = 1:TotalElemCount
    cm = get_cmap(:jet)
    PlotValue=(InputProperty[ElemIdx]-MinValue)/(MaxValue-MinValue)

    face_color = [cm(PlotValue)[1], cm(PlotValue)[2],cm(PlotValue)[3],1.0]

    verts = ((P1[ElemIdx,:],P2[ElemIdx,:],P3[ElemIdx,:]), )
    p3c = PyObject(art3d.Poly3DCollection(verts))
    pycall(ax.add_collection3d, PyAny, p3c)

    # face_color = [0.3, 0.8, 0.3, 0.5]         
    edge_color = [0.2, 0.2, 0.2, 0.0]

    pycall(p3c.set_facecolor, PyAny, face_color)
    pycall(p3c.set_edgecolor, PyAny, edge_color)
    ax.view_init(45, -30)
end

ax.set_aspect("equal")




# SpecificStressStiffness = zeros(TotalElemCount,TotalElemCount)
println("Compiling Stiffness Matrix Function. This may take time if first run")

for ElemIdx = 1:TotalElemCount

    Stress,Strain = TDstressHS(FaultCenter[:,1],FaultCenter[:,2],FaultCenter[:,3],P1[ElemIdx,:],P2[ElemIdx,:],P3[ElemIdx,:],Ss,Ds,Ts,mu,lambda)
    println(ElemIdx)   
    for ElemIdx2 = 1:TotalElemCount
        Stress_i = [Stress[ElemIdx2,1] Stress[ElemIdx2,4] Stress[ElemIdx2,5]
                    Stress[ElemIdx2,4] Stress[ElemIdx2,2] Stress[ElemIdx2,6]
                    Stress[ElemIdx2,5] Stress[ElemIdx2,6] Stress[ElemIdx2,3]]

        TVector = Stress_i * UnitVector_Normal[ElemIdx2,:]
        Stress_Normal = dot(TVector, UnitVector_Normal[ElemIdx2,:])
        Stress_SS = dot(TVector, UnitVector_SS[ElemIdx2,:])
        Stress_Dip = dot(TVector, UnitVector_Dip[ElemIdx2,:])
        Stress_Shear = dot(TVector, UnitVector_Slip[ElemIdx2,:])
        StiffnessMatrix_Normal[ElemIdx2, ElemIdx] = Stress_Normal
        StiffnessMatrix_Shear[ElemIdx2, ElemIdx] = Stress_Shear
        # SpecificStressStiffness[ElemIdx, ElemIdx2] =  Stress[ElemIdx2,6]
    end
end


OutputFileName="Verification/StiffnessCurvedTriangleSSFiner.jld2"
save(OutputFileName, 
"StiffnessMatrixShear", StiffnessMatrix_Shear, "StiffnessMatrixNormal", StiffnessMatrix_Normal, 
"P1",P1,"P2", P2, "P3",P3)
println("Saved File Name: ",OutputFileName)
