using PyPlot
using PyCall
using JLD2
using LinearAlgebra
using Statistics
pygui(true)
include("Functions_Plot.jl")

ResultName="Result"
FileName="Results/" * ResultName * ".jld2"
FileNameInput="Results/" * ResultName * "_Input.jld2"
ResultTime, ResultV, ResultDisp, ResultPressure, Result_NormalStress, History_Theta, =
    load(FileName,"History_Time", "History_V", "History_Disp", "History_Pressure", "History_NormalStress","History_Theta")
ResultV[ResultV.<=0] .= 1e-100



figure(3)
clf()
# PyPlot.plot(ResultTime/60/60/24, log10.(ResultV[:,:]), linewidth=1)
xlabel("Day")
PyPlot.plot(ResultTime/60/60/24/365, log10.(ResultV[:,1:30:end]), linewidth=1)
# #  PyPlot.plot(ResultTime[1:end-1]/60/60/24, log10.(ResultV[1:end-1,1:930]), linewidth=1)
# #  PyPlot.plot(ResultTime[1:end-1]/60/60/24, log10.(ResultV[1:end-1,findall( x -> x == 1, Fault_BulkIndex )]), linewidth=1)


# figure(4)
# clf()
# PyPlot.plot(ResultTime/60/60/24, log10.(History_Theta), linewidth=1)
# xlabel("Day")