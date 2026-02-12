
using PyPlot
using PyCall
using JLD2
using LinearAlgebra
using Statistics
using DelimitedFiles

pygui(true)
include("../../scripts/Functions_Plot.jl")
include("../../CuboidCoupling/Functions_Kuvshinov_Cuboids.jl")

PlotStep = 100
PressurePlotMax = 50e6
# CM_R = get_cmap(:magma_r)   
CM_R = get_cmap(:bwr)    
CM_F = get_cmap(:jet)    


ResultName="Result"
FileName="Results/CuboidResultPlot//" * ResultName * ".jld2"
FileNameInput="Results/CuboidResultPlot//" * ResultName * "_Input.jld2"

Input_Cuboids_File = "Results/CuboidResultPlot/Input_Cuboids.txt"
Input_PorePressure_File = "Results/CuboidResultPlot/Input_PorePressure.txt"
Input_Temperature_File = "Results/CuboidResultPlot/Input_Temperature.txt"


# load Simulation Results
ResultTime, ResultV, ResultDisp, Result_NormalStress, History_Theta, =
load(FileName,"History_Time", "History_V", "History_Disp", "History_NormalStress","History_Theta")
ResultV[ResultV.<=0] .= 1e-100


# Load Cuboids Information
Cuboids_Count, Cuboids_Center, Cuboids_Length = Load_Reservoir_Cuboids(Input_Cuboids_File)
Cuboids_Vertices = Calculate_Cuboids_Vertices(Cuboids_Count, Cuboids_Center, Cuboids_Length)
println("Reservoir cuboid count is:  ", Cuboids_Count)

ExternalStress_TimeArray = readdlm(Input_PorePressure_File, ',' )[1,2:end]
ExternalStress_TimeArray = Float64.(ExternalStress_TimeArray)
TimeArrayCount = length(ExternalStress_TimeArray)

PorePressureChange = readdlm(Input_PorePressure_File, ',')[2:end,2:end]'
PorePressureChange = Float64.(PorePressureChange)

TemperatureChange = readdlm(Input_Temperature_File, ',')[2:end,2:end]'
TemperatureChange = Float64.(TemperatureChange)


# Load Fault Information
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
FaultMass= load(FileNameInput, "FaultMass")

PlotRotation = [59,-83]
Transparent = 1 # 1 for transparent fault plot. 0 for no-transparency
Edge = 0 # 0 for no element boudary. 1 for plotting element boundary
LoadingFaultPlot = 0 # 1 to plot constant velocity faults. 

isdir("3DPlot") || mkdir("3DPlot")
# for PlotStep = 400:400:400
    
    PlotInput=log10.(ResultV[PlotStep,:]); ColorMinMax=[-12, 0] 
    PlotTime = ResultTime[PlotStep]

    PlotIndex = 0
    for Time_idx in eachindex(ExternalStress_TimeArray)
        if PlotTime > ExternalStress_TimeArray[Time_idx]
            PlotIndex = Time_idx 
        end
    end


    Ratio = (PlotTime - ExternalStress_TimeArray[PlotIndex])/(ExternalStress_TimeArray[PlotIndex+1] - ExternalStress_TimeArray[PlotIndex])
 
    P_Gap = PorePressureChange[PlotIndex+1,:] -  PorePressureChange[PlotIndex,:]

    Plot_Pressure = PorePressureChange[PlotIndex,:] + Ratio * P_Gap 
    MaxPressure = PressurePlotMax
    MinPressure = -PressurePlotMax
    Plot_Pressure[Plot_Pressure .> MaxPressure] .=MaxPressure
    Plot_Pressure[Plot_Pressure .< MinPressure] .=MinPressure
    Value = (Plot_Pressure .- MinPressure ) ./ (MaxPressure .- MinPressure)
    Vert1 = Cuboids_Center .- Cuboids_Length /2
    Vert2 = Cuboids_Center .+ Cuboids_Length /2
    figure(1); clf()
    for Cuboididx = 1:Cuboids_Count
        face_color = [CM_R(Value[Cuboididx])[1],CM_R(Value[Cuboididx])[2],CM_R(Value[Cuboididx])[3],1]
        axes = subplot(projection="3d")
        axes.bar3d(
            Vert1[Cuboididx,1], Vert1[Cuboididx,2], -Vert1[Cuboididx,3], Cuboids_Length[Cuboididx,1], Cuboids_Length[Cuboididx,2], -Cuboids_Length[Cuboididx,3],
            color = face_color,
            # alpha=0.5 * Value[Cuboididx] + 0.05,
            alpha=0.3,
            edgecolor=[0.0, 0.0, 0.0,0.3])
    end




    if ColorMinMax == 0 
    MaxValue=maximum(PlotInput)
    MinValue=minimum(PlotInput)
    else
    MaxValue=ColorMinMax[2]
    MinValue=ColorMinMax[1]
    end

    P1, P2, P3 = load(FileNameInput, "P1", "P2", "P3")

    fig = figure(1)
    art3d = PyObject(PyPlot.art3D)
    ax = subplot(projection="3d")
    if Edge == 0 
        edge_color = [0.2, 0.2, 0.2, 0.0]
    else
        edge_color = [0.2, 0.2, 0.2, 0.2]
    end

    for ElemIdx = 1:FaultCount- LoadingFaultCount
        PlotValue=(PlotInput[ElemIdx]-MinValue)/(MaxValue-MinValue)

        if Transparent ==0
            face_color = [CM_F(PlotValue)[1], CM_F(PlotValue)[2],CM_F(PlotValue)[3],1.0]
        else
            face_color = [CM_F(PlotValue)[1], CM_F(PlotValue)[2],CM_F(PlotValue)[3],0.5]
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
    println(PlotTime/60/60/24)
    
    PyPlot.savefig("3DPlot/" * ResultName * "_" * string(PlotStep) * ".png")
# end
