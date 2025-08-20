
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

InputBulkFileName="Input_BulkFaultGeometry.txt"


###### build Principal Stress. Compression Positive. Only Ratio Matters! 
PrincipalStressRatioX = 0.3
PrincipalStressRatioY = 1.0
PrincipalStressRatioZ = 0.6
StressRotationStrike = 30 # degree
StressRotationDip = 0 # degree

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





### Read Bulk Input File

Input_Bulk=readdlm(InputBulkFileName)
Input_Bulk=Input_Bulk[4:end,:]


TotalElemCount = size(Input_Bulk,1)
P1 = Input_Bulk[:,1:3]
P2 = Input_Bulk[:,4:6]
P3 = Input_Bulk[:,7:9]

UnitVector_Normal = zeros(TotalElemCount,3)
UnitVector_Horizontal = zeros(TotalElemCount,3)

for ElemIdx = 1 : TotalElemCount
    P1_i = P1[ElemIdx,:]
    P2_i = P2[ElemIdx,:]
    P3_i = P3[ElemIdx,:]
    UnitVector_Normal_i = cross(P2_i-P1_i, P3_i-P1_i) / norm(cross(P2_i-P1_i, P3_i-P1_i))

    if dot(UnitVector_Normal_i, [1, 0.0001, 0.0001]) >= 0 
    # if angle(UnitVector_Normal_i[1] + UnitVector_Normal_i[2]*im) <= 0 
        P_temp = P1_i
        P1_i = P2_i
        P2_i = P_temp        
    end
    P1[ElemIdx,:] = P1_i
    P2[ElemIdx,:] = P2_i 
    P3[ElemIdx,:] = P3_i 
    UnitVector_Normal[ElemIdx,:] = cross(P2_i-P1_i, P3_i-P1_i) / norm(cross(P2_i-P1_i, P3_i-P1_i))
    UnitVector_Horizontal[ElemIdx,:] = cross(UnitVector_Normal[ElemIdx,:], [0, 0, 1])

end

#     ElemRow = ceil(Int,ElemIdx /TotalElementInRow)
#     OrderInRow = ElemIdx - (ElemRow-1) *TotalElementInRow

#     println("ElemIdx: ", ElemIdx, " ElemRow: ", ElemRow, " OrderInRow: ", OrderInRow)
#     Xvert2 = (OrderInRow - 1) * XElemLength/2
#     Xvert1 = Xvert2 - XElemLength/2
#     Xvert3 = Xvert2 + XElemLength/2
#     Zvert1 = -(ElemRow -1) * ZElemLength
#     Zvert2 = - ElemRow * ZElemLength 
#     if OrderInRow ==1;  Xvert1 = 0.0 ;   end
#     if OrderInRow == TotalElementInRow; Xvert3 = Length;  end
#     if mod(ElemRow,2) == 1
#         if mod(OrderInRow,2) == 1
#             x = [Xvert1 Xvert2 Xvert3]
#             z = [Zvert2 Zvert1 Zvert2]

#         else
#             x = [Xvert1 Xvert2 Xvert3]
#             z = [Zvert1 Zvert2 Zvert1]
#         end    
#     else
#         if mod(OrderInRow,2) == 1
#             x = [Xvert1 Xvert2 Xvert3]
#             z = [Zvert1 Zvert2 Zvert1]

#         else
#             x = [Xvert1 Xvert2 Xvert3]
#             z = [Zvert2 Zvert1 Zvert2]
#         end    
#     end

    
#     # y = x *0  #
#     y = - sind.(x / Length * 180) * Amplitude
#     P1_i = [x[1], y[1], z[1]]
#     P2_i = [x[2], y[2], z[2]]
#     P3_i = [x[3], y[3], z[3]]
#     UnitVector_Normal_i = cross(P2_i-P1_i, P3_i-P1_i) / norm(cross(P2_i-P1_i, P3_i-P1_i))
#     if angle(UnitVector_Normal_i[1] + UnitVector_Normal_i[2]*im) <= 0 
#         P_temp = P1_i
#         P1_i = P2_i
#         P2_i = P_temp        
#     end

#     P1[ElemIdx,:] = P1_i
#     P2[ElemIdx,:] = P2_i
#     P3[ElemIdx,:] = P3_i
#     UnitVector_Normal[ElemIdx,:]= cross(P2_i-P1_i, P3_i-P1_i) / norm(cross(P2_i-P1_i, P3_i-P1_i))
#     # UnitVector_Normal = UnitVector_Normal_i



# end

FaultCenter = (P1+P2+P3)/3
InputProperty = zeros(TotalElemCount)
MaxValue=maximum(InputProperty)
MinValue=minimum(InputProperty)
ArrowLength = 100

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

    ### NormalVector Plot
    ax.quiver(FaultCenter[ElemIdx,1], FaultCenter[ElemIdx,2], FaultCenter[ElemIdx,3], 
        UnitVector_Normal[ElemIdx,1] * ArrowLength, UnitVector_Normal[ElemIdx,2] * ArrowLength,
        UnitVector_Normal[ElemIdx,3] * ArrowLength,
        color="k",arrow_length_ratio=0.2)

        
    ### Horozontal Vector Plot
    ax.quiver(FaultCenter[ElemIdx,1], FaultCenter[ElemIdx,2], FaultCenter[ElemIdx,3], 
        UnitVector_Horizontal[ElemIdx,1] * ArrowLength, UnitVector_Horizontal[ElemIdx,2] * ArrowLength,
        UnitVector_Horizontal[ElemIdx,3] * ArrowLength,
        color="r",arrow_length_ratio=0.2)
end
xlabel("x")
ylabel("y")

ax.set_aspect("equal")



    # ############################# Write Bulk Input #################################
    # ######++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++######
    # #### Bulk File Order 
    # ####  123. P1   456. P2     789. P3     10.Rake
    # ####  11.a      12.b	13.Dc	14.Theta_i	15. V_i     16. Friction_i 17.NormalStress at surface [Pa]  
    # ####  18. NoarmalStress Gradient [Pa] 19. V_Const     20. Minimum Segment Length
    # ######++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++######
    
    # InputFile = [P1  P2  P3 Rake_all a_all b_all Dc_all Thata_i_all V_i_all Friction_i_all NormalStressAtSurface_all NormalStressGradient_all V_Const_all MinimumSegmentLength_all]


    # Rake_all = ones(TotalElemCount) * Rake
    # a_all = ones(TotalElemCount) * a 
    # b_all = ones(TotalElemCount) * b
    # Dc_all = ones(TotalElemCount) * Dc
    # Thata_i_all = ones(TotalElemCount) * Thata_i
    # V_i_all = ones(TotalElemCount) * V_i
    # Friction_i_all = ones(TotalElemCount) * Friction_i
    # NormalStressAtSurface_all = ones(TotalElemCount) * NormalStressAtSurface
    # NormalStressGradient_all = ones(TotalElemCount) * NormalStressGradient
    # V_Const_all = ones(TotalElemCount) * V_Const
     



    # open(OutputFileName, "w") do io
    #     write(io,"SwitchSS/RN\tShearMod\tPoissonRatio\tR_Density\tCrit_TooClose\tTooCloseNormal_Multiplier\tMinimum_NS\n")
    #     writedlm(io,[SwitchSSRN   ShearMod    PoissonRatio     R_Density  Crit_TooClose     TooCloseNormal_Multiplier MinNormalStress])
    #     write(io, "P1X\tP1Y\tP1Z\tP2X\tP2Y\tP2Z\tP3X\tP3Y\tP3Z\tRakeAngle\ta\tb\tDc\tTheta_i\tV_i\tFric_i\tSig0\tSigGrad\tV_Const\tMaxLeng\n")
    #     writedlm(io, InputFile)
    # end;

    # println("Saved File Name: ",OutputFileName)



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
