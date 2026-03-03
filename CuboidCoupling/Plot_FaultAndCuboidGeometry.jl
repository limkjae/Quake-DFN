
using PyPlot
using PyCall
using JLD2
using LinearAlgebra
using Statistics
pygui(true)
include("../scripts/Functions_Plot.jl")
include("Functions_Kuvshinov_Cuboids.jl")


FileNameInput="Input_Discretized.jld2"
Input_Cuboids_File = "CuboidCoupling/Input_Cuboids.txt"

RorT = load(FileNameInput, "RorT")
FaultCount= load(FileNameInput, "FaultCount")
FaultCenter= load(FileNameInput, "FaultCenter")
FaultLengthStrike= load(FileNameInput, "FaultLengthStrike")
FaultLengthDip= load(FileNameInput, "FaultLengthDip")
FaultStrikeAngle= load(FileNameInput, "FaultStrikeAngle")
FaultDipAngle= load(FileNameInput, "FaultDipAngle")
FaultRakeAngle= load(FileNameInput, "FaultRakeAngle")
LoadingFaultCount= load(FileNameInput, "LoadingFaultCount")
Fault_a= load(FileNameInput, "Fault_a")
Fault_b= load(FileNameInput, "Fault_b")
Fault_Dc= load(FileNameInput, "Fault_Dc")
Fault_Theta_i= load(FileNameInput, "Fault_Theta_i")
Fault_V_i= load(FileNameInput, "Fault_V_i")
Fault_Friction_i= load(FileNameInput, "Fault_Friction_i")
Fault_NormalStress= load(FileNameInput, "Fault_NormalStress")
Fault_BulkIndex= load(FileNameInput, "Fault_BulkIndex")

if RorT == "T"
    P1, P2, P3 = load(FileNameInput, "P1", "P2", "P3")
end

Cuboids_Count, Cuboids_Center, Cuboids_Length = Load_Reservoir_Cuboids(Input_Cuboids_File)
Cuboids_Vertices = Calculate_Cuboids_Vertices(Cuboids_Count, Cuboids_Center, Cuboids_Length)
println("Reservoir cuboid count is:  ", Cuboids_Count)


PlotInput = Fault_NormalStress

PlotRotation = [30,-103]
Transparent = 1 # 1 for transparent fault plot. 0 for no-transparency
Edge = 0 # 0 for no element boudary. 1 for plotting element boundary
MinMax_Axis = 0 # 0 for automatically selected axis minimim and maximum 
# MinMax_Axis=[-2000 2000; -2000 2000; -4000 0]
LoadingFaultPlot = 0 # 1 to plot constant velocity faults. 




Vert1 = Cuboids_Center .- Cuboids_Length /2
Vert2 = Cuboids_Center .+ Cuboids_Length /2
figure(1); clf()
for Cuboididx = 1:Cuboids_Count
    face_color = [1, 1, 1, 1]
    axes = subplot(projection="3d")
    axes.bar3d(
        Vert1[Cuboididx,1], Vert1[Cuboididx,2], -Vert1[Cuboididx,3], Cuboids_Length[Cuboididx,1], Cuboids_Length[Cuboididx,2], -Cuboids_Length[Cuboididx,3],
        color = face_color,
        # alpha=0.5 * Value[Cuboididx] + 0.05,
        alpha=0.1,
        edgecolor=[0.0, 0.0, 0.0,0.3])
end

if RorT == "R"
    
    fig= figure(1)

        MaxVaule, MinValue = FaultPlot_3D_Color_General(FaultCenter,FaultLengthStrike, FaultLengthDip,
            FaultStrikeAngle, FaultDipAngle, FaultRakeAngle, PlotInput, 
            PlotRotation, MinMax_Axis, ColorMinMax, Transparent, Edge, LoadingFaultCount)

            ax = subplot(projection="3d")
            PlotTime=ResultTime[PlotStep]/60/60/24
            if ShowDay ==1 
            # ax.text(-5000, 10, 1300, "Day: ",size=10)
            ax.text(DayLocation[1], DayLocation[2], DayLocation[3], PlotTime,size=10)    
            end    
            xlabel("x")
            ylabel("y")
        plotforcbar=  scatter([1,1],[1,1],0.1, [MinValue,MaxVaule], cmap="jet")
        colorbar(plotforcbar, pad=0.15)
        figure(1).canvas.draw()

            ax.set_aspect("equal")

elseif RorT == "T"

    MaxValue=maximum(PlotInput)
    MinValue=minimum(PlotInput)

    fig = figure(1)
    art3d = PyObject(PyPlot.art3D)
    ax = subplot(projection="3d")
    if Edge == 0 
        edge_color = [0.2, 0.2, 0.2, 0.0]
    else
        edge_color = [0.2, 0.2, 0.2, 0.2]
    end

    for ElemIdx = 1:FaultCount- LoadingFaultCount
        cm = get_cmap(:jet)
        PlotValue=(PlotInput[ElemIdx]-MinValue)/(MaxValue-MinValue)

        if Transparent ==0
            face_color = [cm(PlotValue)[1], cm(PlotValue)[2],cm(PlotValue)[3],1.0]
        else
            face_color = [cm(PlotValue)[1], cm(PlotValue)[2],cm(PlotValue)[3],0.5]
        end

        verts = ((P1[ElemIdx,:],P2[ElemIdx,:],P3[ElemIdx,:]), )
        p3c = PyObject(art3d.Poly3DCollection(verts))
        pycall(ax.add_collection3d, PyAny, p3c)
        pycall(p3c.set_facecolor, PyAny, face_color)
        pycall(p3c.set_edgecolor, PyAny, edge_color)
        ax.view_init(PlotRotation[1], PlotRotation[2])
    end
    plotforcbar=  scatter([P1[1,1],P1[2,1]],[P1[1,2],P1[2,2]],0.1, [MinValue,MaxValue], cmap="jet") 
    colorbar(plotforcbar, pad=0.15)
    figure(1).canvas.draw()
    ax.set_aspect("equal")
end
