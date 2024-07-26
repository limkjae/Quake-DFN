using PyPlot
using PyCall
using JLD2
using LinearAlgebra
using Statistics

pygui(true)
include("Functions_Plot.jl")

PlotBulkIndex = 2
TimeIndex = 800

ResultName="Result"
FileName="Results/" * ResultName * ".jld2"
FileNameInput="Results/" * ResultName * "_Input.jld2"

ResultTime, ResultV, ResultDisp, ResultPressure, Result_NormalStress, Result_Theta =
    load(FileName,"History_Time", "History_V", "History_Disp", "History_Pressure", "History_NormalStress", "History_Theta")
# ResultV[ResultV.<=1e-12] .= 1e-12

# figure(10); PyPlot.plot(log10.(ResultV[:,1:30:end])); xlabel("Record Step")

FaultCount= load(FileNameInput, "FaultCount")
FaultCenter= load(FileNameInput, "FaultCenter")
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
ShearModulus= load(FileNameInput,"ShearModulus")

TheBulkIndexes = findall(x->x==PlotBulkIndex, Fault_BulkIndex)
MinX = minimum(FaultCenter[TheBulkIndexes,1])
MinY = minimum(FaultCenter[TheBulkIndexes,2])
SecondMinX = maximum(FaultCenter[TheBulkIndexes,1])
SecondMinY = maximum(FaultCenter[TheBulkIndexes,2])
for ElementIndex in TheBulkIndexes
    if FaultCenter[ElementIndex,1] > MinX && FaultCenter[ElementIndex,1] < SecondMinX
        SecondMinX = FaultCenter[ElementIndex,1]
    end
    if FaultCenter[ElementIndex,2] > MinY && FaultCenter[ElementIndex,2] < SecondMinY
        SecondMinY = FaultCenter[ElementIndex,2]
    end
end

GapX = SecondMinX - MinX
Xindex = round.(Int, (FaultCenter[TheBulkIndexes,1] .- MinX) / GapX .+ 1 )

GapY = SecondMinY - MinY
Yindex = round.(Int, (FaultCenter[TheBulkIndexes,2] .- MinY) / GapY .+ 1 )


ResultV[ResultV.<=1e-15] .= 1e-15

PlotTime=ResultTime[TimeIndex]/60/60/24
println(PlotTime)

ContourX = zeros(maximum(Xindex),maximum(Yindex))
ContourY = zeros(maximum(Xindex),maximum(Yindex))
ContourPlot= zeros(maximum(Xindex),maximum(Yindex))
levels = []
for i in eachindex(TheBulkIndexes)
    ContourX[Xindex[i],Yindex[i]] = FaultCenter[TheBulkIndexes[i],1]    
    ContourY[Xindex[i],Yindex[i]] = FaultCenter[TheBulkIndexes[i],2]
    ContourPlot[Xindex[i],Yindex[i]] = log10(ResultV[TimeIndex,TheBulkIndexes[i]]); levels = collect(-15:0.1:1)
    # ContourPlot[Xindex[i],Yindex[i]] = ResultDisp[TimeIndex,i]; levels = collect(0:0.01:1)
end
# figure(2)


fig1 = contourf(ContourX,ContourY,ContourPlot, levels, cmap="jet")
plt.pause(0.02)
plt.show()
