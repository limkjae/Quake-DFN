
using PyPlot
using PyCall
using JLD2
using LinearAlgebra
using DelimitedFiles
pygui(true)

ResultName="Result"


MarkerSizeAdjustment = -2
MergeTimeCriteria = 10.0
MergeDistanceCriteria = 5000.0
PlotRotation=[25,-30]
DispRateCrits = 1e-2



FileName="Results/" * ResultName * ".jld2"
FileNameInput="Results/" * ResultName * "_Input.jld2"

include("Functions_Plot.jl")


EventTime_Fragment, EventLocation_Fragment, EventMoment_Fragment, EventCount_Fragment = 
     get_event_fragments(FileName, FileNameInput)

EventTime_Merged, EventMomentMagnitude,  EventLocation_Merged =  
     merge_fragments(EventTime_Fragment, EventLocation_Fragment, EventMoment_Fragment, EventCount_Fragment, MergeTimeCriteria, MergeDistanceCriteria)




######################## Plots ##########################

figure(1)
scatter(EventTime_Merged/60/60/24, EventMomentMagnitude,s= 6 .^ (EventMomentMagnitude .+ MarkerSizeAdjustment), facecolors="none", edgecolor="b")
ylim([2,8])


figure(2)
scatter3D(EventLocation_Merged[:,1], EventLocation_Merged[:,2],-EventLocation_Merged[:,3],s= 6 .^ (EventMomentMagnitude .+ MarkerSizeAdjustment), facecolors="none", edgecolor="b")

figure(3)
scatter(EventLocation_Merged[:,1], EventLocation_Merged[:,2],s= 6 .^ (EventMomentMagnitude .+ MarkerSizeAdjustment), facecolors="none", edgecolor="b")



##################### Save Catalog #######################
OutputFileName_MomentVsTime="Results/Catalog_" *ResultName
OutputMomentVsTime = [EventTime_Merged'; EventMomentMagnitude'; EventLocation_Merged']'

open(OutputFileName_MomentVsTime, "w") do io
     writedlm(io, OutputMomentVsTime)
     end
          






     

