
using PyPlot
using PyCall
using JLD2
pygui(true)

ResultName="Result"

include("Functions_Plot.jl")
#include("Results")
FileName="Results/" * ResultName * ".jld2"
FileNameInput="Results/" * ResultName * "_Input.jld2"


FaultCount= load(FileNameInput, "FaultCount")
FaultCenter= load(FileNameInput, "FaultCenter")
FaultLengthStrike= load(FileNameInput, "FaultLengthStrike")
FaultLengthDip= load(FileNameInput, "FaultLengthDip")
FaultStrikeAngle= load(FileNameInput, "FaultStrikeAngle")
FaultDipAngle= load(FileNameInput, "FaultDipAngle")
FaultLLRR= load(FileNameInput, "FaultLLRR")
LoadingFaultCount= load(FileNameInput, "LoadingFaultCount")
Fault_a= load(FileNameInput, "Fault_a")
Fault_b= load(FileNameInput, "Fault_b")
Fault_Dc= load(FileNameInput, "Fault_Dc")
Fault_Theta_i= load(FileNameInput, "Fault_Theta_i")
Fault_V_i= load(FileNameInput, "Fault_V_i")
Fault_Friction_i= load(FileNameInput, "Fault_Friction_i")
Fault_NormalStress= load(FileNameInput, "Fault_NormalStress")


# PlotInput =Fault_a - Fault_b
PlotInput =Fault_Dc


PlotRotation=[45,-45]
Transparent=1
Edge= 0 
MinMax_Axis=0 # automatically detect max and min 
# MinMax_Axis=[-3000 3000; -3000 3000; -4000 0]

ColorMinMax=0
# ColorMinMax=[-10,-1]

LoadingFaultPlot=0 # 1 to plot constant velocity faults. 





if LoadingFaultPlot==1
    LoadingFaultCount=0
end
    

figure(1)
clf()
MaxVaule, MinValue = FaultPlot_3D_Color_General(FaultCenter,FaultLengthStrike, FaultLengthDip,
    FaultStrikeAngle, FaultDipAngle, FaultLLRR, PlotInput, 
    PlotRotation, MinMax_Axis, ColorMinMax, Transparent, Edge, LoadingFaultCount)

# figure(1)
plotforcbar=  scatter([1,1],[1,1],0.1, [MinValue,MaxVaule], cmap="jet")
colorbar(plotforcbar, pad=0.15)
figure(1).canvas.draw()

    
