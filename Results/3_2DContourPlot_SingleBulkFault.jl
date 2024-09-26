using PyPlot
using PyCall
using JLD2
using LinearAlgebra
using Statistics

pygui(true)
include("Functions_Plot.jl")

PlotBulkIndex = 2
TimeIndex = 100

XPlotIndex = 2
YPlotIndex = 3
    
for OneLoop = 1:1
    println(OneLoop)
    ResultName="Result"
    FileName="Results/" * ResultName * ".jld2"
    FileNameInput="Results/" * ResultName * "_Input.jld2"

    ResultTime, ResultV, ResultDisp, Result_NormalStress, Result_Theta =
        load(FileName,"History_Time", "History_V", "History_Disp", "History_NormalStress", "History_Theta")

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
    MinX = minimum(FaultCenter[TheBulkIndexes,XPlotIndex])
    MinY = minimum(FaultCenter[TheBulkIndexes,YPlotIndex])
    SecondMinX = maximum(FaultCenter[TheBulkIndexes,XPlotIndex])
    SecondMinY = maximum(FaultCenter[TheBulkIndexes,YPlotIndex])
    for ElementIndex in TheBulkIndexes
        if FaultCenter[ElementIndex,XPlotIndex] > MinX && FaultCenter[ElementIndex,XPlotIndex] < SecondMinX
            SecondMinX = FaultCenter[ElementIndex,XPlotIndex]
        end
        if FaultCenter[ElementIndex,YPlotIndex] > MinY && FaultCenter[ElementIndex,YPlotIndex] < SecondMinY
            SecondMinY = FaultCenter[ElementIndex,YPlotIndex]
        end
    end

    GapX = SecondMinX - MinX
    if GapX == 0
        println("There is no distribution in PlotX direction. Consider change XPlotIndex")
        break
    end
    Xindex = round.(Int, (FaultCenter[TheBulkIndexes,XPlotIndex] .- MinX) / GapX .+ 1 )

    GapY = SecondMinY - MinY
    if GapY == 0
        println("There is no distribution in PlotY direction. Consider change YPlotIndex")
        break 
    end
    Yindex = round.(Int, (FaultCenter[TheBulkIndexes,YPlotIndex] .- MinY) / GapY .+ 1 )


    ResultV[ResultV.<=1e-15] .= 1e-15

    PlotTime=ResultTime[TimeIndex]/60/60/24
    println("Plot time is: ", PlotTime)

    ContourX = zeros(maximum(Xindex),maximum(Yindex))
    ContourY = zeros(maximum(Xindex),maximum(Yindex))
    ContourPlot= zeros(maximum(Xindex),maximum(Yindex))
    levels = []
    for i in eachindex(TheBulkIndexes)
        ContourX[Xindex[i],Yindex[i]] = FaultCenter[TheBulkIndexes[i],XPlotIndex]    
        ContourY[Xindex[i],Yindex[i]] = FaultCenter[TheBulkIndexes[i],YPlotIndex]
        ContourPlot[Xindex[i],Yindex[i]] = log10(ResultV[TimeIndex,TheBulkIndexes[i]]); levels = collect(-15:0.1:1)
        # ContourPlot[Xindex[i],Yindex[i]] = ResultDisp[TimeIndex,i]; levels = collect(0:0.01:1)
    end
    # figure(2)


    fig1 = contourf(ContourX,ContourY,ContourPlot, levels, cmap="jet")
    plt.pause(0.02)
    plt.show()
end

