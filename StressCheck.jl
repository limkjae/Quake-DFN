

using PyPlot
using PyCall
using DelimitedFiles
using JLD2
using LinearAlgebra
using Printf
using SpecialFunctions
using StaticArrays
pygui(true)

include("Functions_Solvers.jl")
include("Functions_RSFDFN3DMain_D.jl")
include("Results/Functions_Plot.jl")
include("QuickParameterAdjust.jl")

SlipBulkIndex = 1

SaveResultFileName="Results/Result.jld2"
SaveInputInfoFileName="Results/Result_Input.jld2" 

LoadingInputFileName="Input_Discretized.jld2" 
StiffnessMatrixShear= load(LoadingInputFileName, "StiffnessMatrixShear")
StiffnessMatrixNormal= load(LoadingInputFileName, "StiffnessMatrixNormal")

FaultCenter= load(LoadingInputFileName, "FaultCenter")
ShearModulus= load(LoadingInputFileName, "ShearModulus")
FaultLengthStrike= load(LoadingInputFileName, "FaultLengthStrike")
FaultLengthDip= load(LoadingInputFileName, "FaultLengthDip")
FaultStrikeAngle= load(LoadingInputFileName, "FaultStrikeAngle")
FaultDipAngle= load(LoadingInputFileName, "FaultDipAngle")
FaultLLRR= load(LoadingInputFileName, "FaultLLRR")


Fault_BulkIndex = Int.(load(LoadingInputFileName, "Fault_BulkIndex"))

ShearSlip = zeros(length(Fault_BulkIndex))
for ElemIdx in eachindex(Fault_BulkIndex)
    if Fault_BulkIndex[ElemIdx] == SlipBulkIndex
        ShearSlip[ElemIdx] = 1
    end
end

ShearStressChange = StiffnessMatrixShear * ShearSlip
NormalStressChange = StiffnessMatrixNormal * ShearSlip
for ElemIdx in eachindex(Fault_BulkIndex)
    if Fault_BulkIndex[ElemIdx] == SlipBulkIndex
        ShearStressChange[ElemIdx] = 0.0
        NormalStressChange[ElemIdx] = 0.0
    end
end





# PlotInput = ShearStressChange; ColorMinMax=[-1e6,1e6]
PlotInput = NormalStressChange;ColorMinMax=[-1e6,1e6]

PlotRotation = [35,-30]
Transparent = 0 # 1 for transparent fault plot. 0 for no-transparency
Edge = 0 # 0 for no element boudary. 1 for plotting element boundary
MinMax_Axis = 0 # 0 for automatically selected axis minimim and maximum 
# MinMax_Axis=[-2000 2000; -2000 2000; -4000 0]
LoadingFaultPlot = 0 # 1 to plot constant velocity faults. 

ShowDay = 0 # If 1, day is shown in the location 
DayLocation = [0,0,1000]
LoadingFaultCount = 0

figure(1)
clf()
MaxVaule, MinValue = FaultPlot_3D_Color_General(FaultCenter,FaultLengthStrike, FaultLengthDip,
    FaultStrikeAngle, FaultDipAngle, FaultLLRR, PlotInput, 
    PlotRotation, MinMax_Axis, ColorMinMax, Transparent, Edge, LoadingFaultCount)

    ax = subplot(projection="3d")
    xlabel("x")
    ylabel("y")
plotforcbar=  scatter([1,1],[1,1],0.1, [MinValue,MaxVaule], cmap="jet")
colorbar(plotforcbar, pad=0.15)
figure(1).canvas.draw()

