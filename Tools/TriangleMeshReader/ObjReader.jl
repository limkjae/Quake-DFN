
using DelimitedFiles
using PyPlot
using PyCall
using LinearAlgebra
using Statistics
using JLD2

pygui(true)

FileName = "tripleJunction.obj"


InputOBJFileName="Tools/TriangleMeshReader/$(FileName)"
InputBulkFileName="Input_BulkFaultGeometry.txt"


SwitchSSRN = 0
ShearMod = 20e9
PoissonRatio = 0.25
R_Density=2670.0
Crit_TooClose= 1.05
TooCloseNormal_Multiplier = 0.6
MinNormalStress = 2e6


Rake = 0
a = 0.005
b = 0.01
Dc = 3e-4
Thata_i = 1e5
V_i = 1e-15
Friction_i = 0.6
NormalStressAtSurface = 10e6
NormalStressGradient = 0
V_Const = 0.0
MinimumSegmentLength = 100000.0






MeshOBJ_Raw=readdlm(InputOBJFileName)
VertexCount = 0
FaceCount = 0
VertexNormalCount = 0
Vertex = zeros(1,3)
VertexNormal = zeros(1,3)
Face = zeros(Int, 1,3)
for RowI = 1:length(MeshOBJ_Raw[:,1])
    if MeshOBJ_Raw[RowI,1] == "v"
        
        Vertex = [Vertex;    vec(MeshOBJ_Raw[RowI,2:4])']

    elseif MeshOBJ_Raw[RowI,1] == "vn"

        VertexNormal= [Vertex;    vec(MeshOBJ_Raw[RowI,2:4])']

    elseif MeshOBJ_Raw[RowI,1] == "f"

        Face = [Face; [parse(Int64,split(MeshOBJ_Raw[RowI,2], '/')[1]), parse(Int64,split(MeshOBJ_Raw[RowI,3], '/')[1]), parse(Int64,split(MeshOBJ_Raw[RowI,4], '/')[1])]']

    end

end

Vertex = Vertex[2:end,:] * 100
VertexNormal = VertexNormal[2:end,:]
Face = Face[2:end,:]
Face[2,1]
TotalElemCount = length(Face[:,1])

InputProperty = rand(TotalElemCount)
MaxValue=maximum(InputProperty)
MinValue=minimum(InputProperty)

figure(1)
fig = figure(1)
clf()
art3d = PyObject(PyPlot.art3D)
ax = subplot(projection="3d")
P1 = zeros(TotalElemCount,3)
P2 = zeros(TotalElemCount,3)
P3 = zeros(TotalElemCount,3)
for ElemIdx = 1:TotalElemCount
    cm = get_cmap(:jet)
    PlotValue=(InputProperty[ElemIdx]-MinValue)/(MaxValue-MinValue)

    face_color = [cm(PlotValue)[1], cm(PlotValue)[2],cm(PlotValue)[3],0.0]
    VertElem = [Face[ElemIdx,1], Face[ElemIdx,2], Face[ElemIdx,3]]
    verts = ((Vertex[VertElem[1],:],Vertex[VertElem[2],:],Vertex[VertElem[3],:]), )
    # verts = ((P1[ElemIdx,:],P2[ElemIdx,:],P3[ElemIdx,:]), )
    p3c = PyObject(art3d.Poly3DCollection(verts))
    pycall(ax.add_collection3d, PyAny, p3c)

    # face_color = [0.3, 0.8, 0.3, 0.5]         
    edge_color = [0.2, 0.2, 0.2, 0.5]

    pycall(p3c.set_facecolor, PyAny, face_color)
    pycall(p3c.set_edgecolor, PyAny, edge_color)
    ax.view_init(45, -30)
    P1[ElemIdx,:] = Vertex[VertElem[1],:]
    P2[ElemIdx,:] = Vertex[VertElem[2],:]
    P3[ElemIdx,:] = Vertex[VertElem[3],:]
end
ax.set_aspect("equal")

Input_Bulk = zeros(TotalElemCount,20)
Input_Bulk[:,1:3] = P1
Input_Bulk[:,4:6] = P2
Input_Bulk[:,7:9] = P3
Input_Bulk[:,10] .= 0
Input_Bulk[:,11] .= a
Input_Bulk[:,12] .= b
Input_Bulk[:,13] .= Dc
Input_Bulk[:,14] .= Thata_i
Input_Bulk[:,15] .= V_i
Input_Bulk[:,16] .= Friction_i
Input_Bulk[:,17] .= NormalStressAtSurface
Input_Bulk[:,18] .= NormalStressGradient
Input_Bulk[:,19] .= V_Const
Input_Bulk[:,20] .= MinimumSegmentLength


Input_BulkHeader = [SwitchSSRN ShearMod PoissonRatio R_Density Crit_TooClose TooCloseNormal_Multiplier MinNormalStress]

    ############################# Write Bulk Input #################################
    ######++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++######
    #### Bulk File Order 
    ####  123. P1   456. P2     789. P3     10.Rake
    ####  11.a      12.b	13.Dc	14.Theta_i	15. V_i     16. Friction_i 17.NormalStress at surface [Pa]  
    ####  18. NoarmalStress Gradient [Pa] 19. V_Const     20. Minimum Segment Length
    ######++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++######

    open(InputBulkFileName, "w") do io
        write(io,"SwitchSS/RN\tShearMod\tPoissonRatio\tR_Density\tCrit_TooClose\tTooCloseNormal_Multiplier\tMinimum_NS\n")
        writedlm(io, Input_BulkHeader)
        write(io, "P1_x\tP1_y\tP1_z\tP2_x\tP2_y\tP2_z\tP3_x\tP3_y\tP3_z\tRake\ta\tb\tDc\tTheta_i\tV_i\tFric_i\tSig0\tSigGrad\tV_Const\tMaxLeng\n")
        writedlm(io, Input_Bulk)
    end

area = norm.(cross.(eachrow(Input_Bulk[:,1:3] - Input_Bulk[:,4:6]), eachrow(Input_Bulk[:,1:3] - Input_Bulk[:,7:9])))/2
AverageArea = mean(area)
TooSmall = findall(area .< AverageArea /100)
if length(TooSmall)>1
    println(TooSmall)
    println("The area of above elements seems to be too small, please check mesh")
    
end