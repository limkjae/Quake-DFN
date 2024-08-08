
using PyPlot
using PyCall
using JLD2
using LinearAlgebra
using DelimitedFiles
pygui(true)

ResultName="Result"

InjectionLocation= [0.0,0.0,5000]
MarkerSizeAdjustment = 0
MergeTimeCriteria = 1.0
MergeDistanceCriteria = 2000.0
PlotRotation=[25,-30]
DispRateCrits = 1e-2



FileName="Results/" * ResultName * ".jld2"
FileNameInput="Results/" * ResultName * "_Input.jld2"

include("Functions_Plot.jl")

ResultTime=load(FileName,"History_Time")
ResultDisp=load(FileName,"History_Disp")
ResultV=load(FileName,"History_V")
FaultCenter =load(FileNameInput,"FaultCenter")
FaultCount =load(FileNameInput,"FaultCount")
FaultLengthStrike =load(FileNameInput, "FaultLengthStrike")
FaultLengthDip =load(FileNameInput, "FaultLengthDip")
Fault_BulkIndex= load(FileNameInput, "Fault_BulkIndex")


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
ShearModulus=load(FileNameInput, "ShearModulus")


FragmentMomentAccum = zeros(size(ResultDisp,1),size(ResultDisp,2))
FragmentMomentRate = zeros(size(ResultDisp,1),size(ResultDisp,2))
DispRate = zeros(size(ResultDisp,1),size(ResultDisp,2))

EventCount_Fragment=1
EventTime_Fragment=[0;]
EventMoment_Fragment=[0;]
EventLocation_Fragment=[0 0 0]


log10.(ResultV[1:end-1,findall( x -> x == 1, Fault_BulkIndex )])

for FaultIdx in eachindex(ResultDisp[1,:])
     EventOnOrOff=0
     for TimeIdx = 2:length(ResultDisp[:,1])

          FragmentMomentAccum[TimeIdx,FaultIdx]=ResultDisp[TimeIdx,FaultIdx] * 
               FaultLengthStrike[FaultIdx] * FaultLengthDip[FaultIdx]  * ShearModulus           
          DispRate[TimeIdx,FaultIdx]= (ResultDisp[TimeIdx,FaultIdx] - ResultDisp[TimeIdx-1,FaultIdx] ) / (ResultTime[TimeIdx] - ResultTime[TimeIdx-1])

          if DispRate[TimeIdx, FaultIdx] > DispRateCrits
               if EventOnOrOff == 0
                    EventOnOrOff=1
                    EventCount_Fragment=EventCount_Fragment+1
                    EventTime_Fragment=[EventTime_Fragment; ResultTime[TimeIdx]]
                    EventLocation_Fragment=[EventLocation_Fragment; FaultCenter[FaultIdx,:]']
                    EventMoment_Fragment=[EventMoment_Fragment; FragmentMomentAccum[TimeIdx,FaultIdx] - FragmentMomentAccum[TimeIdx-1,FaultIdx]]
               else                    
                    EventMoment_Fragment[EventCount_Fragment] = EventMoment_Fragment[EventCount_Fragment] + FragmentMomentAccum[TimeIdx,FaultIdx] - FragmentMomentAccum[TimeIdx-1,FaultIdx]
               end
          else
               EventOnOrOff = 0
          end
     end
end


EventCount_Fragment=EventCount_Fragment-1
EventTime_Fragment=EventTime_Fragment[2:end]
EventLocation_Fragment=EventLocation_Fragment[2:end,:]
EventMoment_Fragment=EventMoment_Fragment[2:end]




TotalMomentAccum=sum(FragmentMomentAccum,dims=2)

######### Rearrange the fragments by time order ##########
SortOrder=sortperm(EventTime_Fragment)
EventTime_Fragment=[EventTime_Fragment[SortOrder];]
EventLocation_Fragment=EventLocation_Fragment[SortOrder , :]
EventMoment_Fragment=[EventMoment_Fragment[SortOrder];]



EventRelLocFromSource=zeros(size(EventLocation_Fragment))
for i in eachindex(EventRelLocFromSource[:,1])
EventRelLocFromSource[i,:]=EventLocation_Fragment[i,:] - InjectionLocation
end

EventDistance=sqrt.(sum(EventRelLocFromSource.^2 , dims=2))



SegEventBackIdx=1
EventCount_Merged=0
EventTime_Merged=[0;]
EventLocation_Merged=[0 0 0]

EventNumberMerged_forFragment=zeros(Int, EventCount_Fragment)
SumMoment=0.0
CurrentBulkNumber=0
for EventFragmentIdx=1:EventCount_Fragment

     Terminate = 0 
     SegEventBackIdx = 1
     while Terminate == 0 
          CheckFragmentIdx = EventFragmentIdx - SegEventBackIdx # check fragment until given time 
          if EventFragmentIdx == 1 # initial setting
               EventCount_Merged = 1
               EventTime_Merged = EventTime_Fragment[EventFragmentIdx]
               EventLocation_Merged = EventLocation_Fragment[EventFragmentIdx, :]'
               EventNumberMerged_forFragment[EventFragmentIdx] = EventCount_Merged
               Terminate = 1 
          elseif  EventTime_Fragment[EventFragmentIdx] -  EventTime_Fragment[CheckFragmentIdx] <  MergeTimeCriteria
               # if within the time range
               if norm(EventLocation_Fragment[EventFragmentIdx,:] - EventLocation_Fragment[CheckFragmentIdx,:]) < MergeDistanceCriteria
                    # if within the distance range -> This is same event to this event (time and distance satisfied)
                    EventNumberMerged_forFragment[EventFragmentIdx] = EventNumberMerged_forFragment[CheckFragmentIdx] 
                    Terminate = 1 
               else
                    #this event is within the time range but not in distance range. keep checking
                    SegEventBackIdx = SegEventBackIdx + 1
               end
          else
               # Now time range exceeded without finding merged event
               # This is a new event
               EventCount_Merged = EventCount_Merged + 1
               EventTime_Merged = [EventTime_Merged; EventTime_Fragment[EventFragmentIdx]]
               EventLocation_Merged = [EventLocation_Merged ; EventLocation_Fragment[EventFragmentIdx, :]']
               EventNumberMerged_forFragment[EventFragmentIdx] = EventCount_Merged
               Terminate = 1 
          end
     end
end


# EventTime_Merged = EventTime_Merged[2:end]
# EventLocation_Merged = EventLocation_Merged[2:end,:]


println("new event count = ", EventCount_Merged)

EventMoment_Bulk = zeros(EventCount_Merged)

for EventFragmentIdx=1:EventCount_Fragment
     MergedIdxofThisFragment = EventNumberMerged_forFragment[EventFragmentIdx]
     EventMoment_Bulk[MergedIdxofThisFragment] = EventMoment_Bulk[MergedIdxofThisFragment] + EventMoment_Fragment[EventFragmentIdx]
end
EventMomentMagnitude = (log10.(EventMoment_Bulk) .- 9.05) ./ 1.5


EventDistanceFromSource_Bulk = zeros(EventCount_Merged)
for EventFragmentIdx=1:EventCount_Merged
     EventDistanceFromSource_Bulk[EventFragmentIdx] = norm(EventLocation_Merged[EventFragmentIdx,:] - InjectionLocation)

end


# MarkerSizeAdjustment = 0

# figure(2) 
# clf()
# # scatter(log10.(EventTime_Merged/60/60/24), EventDistanceFromSource_Bulk, s= 6 .^ (EventMomentMagnitude .+ MarkerSizeAdjustment), facecolors="none", edgecolor="k" )

# scatter(EventTime_Merged/60/60/24, EventDistanceFromSource_Bulk, s= 6 .^ (EventMomentMagnitude .+ MarkerSizeAdjustment), facecolors="none", edgecolor="k" )


# figure(3)
# clf()
# scatter(EventTime_Merged/60/60/24, EventLocation_Merged[:,2], s= 6 .^ (EventMomentMagnitude .+ MarkerSizeAdjustment), facecolors="none", edgecolor="k" )

# figure(4)
# clf()
# scatter(EventTime_Merged/60/60/24, -EventLocation_Merged[:,3], s= 6 .^ (EventMomentMagnitude .+ MarkerSizeAdjustment), facecolors="none", edgecolor="k" )

figure(5)
# clf()
# scatter(log10.(EventTime_Merged/60/60/24), EventMomentMagnitude, s= 6 .^ (EventMomentMagnitude .+ MarkerSizeAdjustment), facecolors="none", edgecolor="k" )
scatter(EventTime_Merged/60/60/24, EventMomentMagnitude,s= 6 .^ (EventMomentMagnitude .+ MarkerSizeAdjustment), facecolors="none", edgecolor="b")
ylim([2,8])


figure(6)
# clf()
# scatter(log10.(EventTime_Merged/60/60/24), EventMomentMagnitude, s= 6 .^ (EventMomentMagnitude .+ MarkerSizeAdjustment), facecolors="none", edgecolor="k" )
scatter3D(EventLocation_Merged[:,1], EventLocation_Merged[:,2],-EventLocation_Merged[:,3],s= 6 .^ (EventMomentMagnitude .+ MarkerSizeAdjustment), facecolors="none", edgecolor="b")

figure(7)
# clf()
# scatter(log10.(EventTime_Merged/60/60/24), EventMomentMagnitude, s= 6 .^ (EventMomentMagnitude .+ MarkerSizeAdjustment), facecolors="none", edgecolor="k" )
scatter(EventLocation_Merged[:,1], EventLocation_Merged[:,2],s= 6 .^ (EventMomentMagnitude .+ MarkerSizeAdjustment), facecolors="none", edgecolor="b")
# xlim([-400,400])



OutputFileName_MomentVsTime="Results/Catalog_" *ResultName
OutputMomentVsTime = [EventTime_Merged'; EventMomentMagnitude'; EventLocation_Merged'; EventDistanceFromSource_Bulk']'

open(OutputFileName_MomentVsTime, "w") do io
     writedlm(io, OutputMomentVsTime)
     end
          






     

