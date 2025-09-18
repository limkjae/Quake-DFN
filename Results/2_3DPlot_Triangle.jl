
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

ResultTime, ResultV, ResultDisp, Result_NormalStress, History_Theta, =
load(FileName,"History_Time", "History_V", "History_Disp", "History_NormalStress","History_Theta")
ResultV[ResultV.<=0] .= 1e-100

#######################################################################################
############################### Figure Configuration ##################################

# figure(10); clf(); PyPlot.plot(log10.(ResultV[:,1:50:end])); xlabel("Record Step")
PlotStep = 450

PlotRotation = [43,-44]
Transparent = 0 # 1 for transparent fault plot. 0 for no-transparency
Edge = 0 # 0 for no element boudary. 1 for plotting element boundary
MinMax_Axis = 0 # 0 for automatically selected axis minimim and maximum 
# MinMax_Axis=[-2000 2000; -2000 2000; -4000 0]
LoadingFaultPlot = 0 # 1 to plot constant velocity faults. 

ShowDay = 0 # If 1, day is shown in the location 
DayLocation = [0,0,1000]
# MaxVLog = log10(maximum(ResultV[:,:]))
################## What to Plot ? ###################
PlotInput=log10.(ResultV[PlotStep,:]); ColorMinMax=[-12, 0] 
# PlotInput=log10.(ResultV[PlotStep,:]); ColorMinMax=[MaxVLog-3,MaxVLog] 
# PlotInput= Result_NormalStress[PlotStep,1:end-4] ; ColorMinMax=0
# PlotInput=log10.(ResultPressure[PlotStep,:]); ColorMinMax=[3,6]
# PlotInput= Result_NormalStress[PlotStep,:] -  Fault_NormalStress; ColorMinMax=[-1e6,1e6]
# PlotInput=ResultDisp[PlotStep,:]; ColorMinMax= [0.03, 0.08] 

#############---------------------------#############



##############+++++++++++++++++++++++++++++++++++++++++++++++++++++++++################
#######################################################################################



P1, P2, P3 = load(FileNameInput, "P1", "P2", "P3")
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
TotalElemCount = FaultCount
if ColorMinMax == 0 
MaxValue=maximum(PlotInput)
MinValue=minimum(PlotInput)
else
MaxValue=ColorMinMax[2]
MinValue=ColorMinMax[1]
end
figure(1)
fig = figure(1)
clf()
art3d = PyObject(PyPlot.art3D)
ax = subplot(projection="3d")
for ElemIdx = 1:TotalElemCount- LoadingFaultCount
    cm = get_cmap(:jet)
    PlotValue=(PlotInput[ElemIdx]-MinValue)/(MaxValue-MinValue)

    face_color = [cm(PlotValue)[1], cm(PlotValue)[2],cm(PlotValue)[3],1.0]

    verts = ((P1[ElemIdx,:],P2[ElemIdx,:],P3[ElemIdx,:]), )
    p3c = PyObject(art3d.Poly3DCollection(verts))
    pycall(ax.add_collection3d, PyAny, p3c)

    # face_color = [0.3, 0.8, 0.3, 0.5]         
    edge_color = [0.2, 0.2, 0.2, 0.0]

    pycall(p3c.set_facecolor, PyAny, face_color)
    pycall(p3c.set_edgecolor, PyAny, edge_color)
    ax.view_init(PlotRotation[1], PlotRotation[2])
end

ax.set_aspect("equal")





# if LoadingFaultPlot==1
# LoadingFaultCount=0
# end

# figure(1)
# clf()
# MaxVaule, MinValue = FaultPlot_3D_Color_General(FaultCenter,FaultLengthStrike, FaultLengthDip,
#     FaultStrikeAngle, FaultDipAngle, FaultRakeAngle, PlotInput, 
#     PlotRotation, MinMax_Axis, ColorMinMax, Transparent, Edge, LoadingFaultCount)

#     ax = subplot(projection="3d")
#     PlotTime=ResultTime[PlotStep]/60/60/24
#     if ShowDay ==1 
#     # ax.text(-5000, 10, 1300, "Day: ",size=10)
#     ax.text(DayLocation[1], DayLocation[2], DayLocation[3], PlotTime,size=10)    
#     end    
#     xlabel("x")
#     ylabel("y")
# plotforcbar=  scatter([1,1],[1,1],0.1, [MinValue,MaxVaule], cmap="jet")
# colorbar(plotforcbar, pad=0.15)
# figure(1).canvas.draw()




# if Animation_Save == 1

#     # PlotStepI=PlotStep
#     isdir("3DPlot") || mkdir("3DPlot")
#     for i=StepBegin:StepInterval:StepEnd
#         PlotStep = i

#         clf()
#         PlotInput=log10.(ResultV[PlotStep,:])

#         MaxVaule, MinValue = FaultPlot_3D_Color_General(FaultCenter,FaultLengthStrike, FaultLengthDip,
#         FaultStrikeAngle, FaultDipAngle, FaultRakeAngle, PlotInput, 
#         PlotRotation, MinMax_Axis, ColorMinMax, Transparent, Edge, LoadingFaultCount)
#         figure(1).canvas.draw()
#         ax = subplot(projection="3d")
#         PlotTime=ResultTime[PlotStep]/60/60/24
#         if ShowDay ==1 
#         # ax.text(-5000, 10, 1300, "Day: ",size=10)
#         ax.text(DayLocation[1], DayLocation[2], DayLocation[3], PlotTime,size=10)    
#         end    
#         figure(1).canvas.draw()

#         println(PlotStep)
        
#         PyPlot.savefig("3DPlot/" * ResultName * "_" * string(i) * ".png")
#     end

# end



# # #Single element locator
# # figure(1)
# SelectedElements = [1275,1277,1278,1274]
# # clf()
# MaxVaule, MinValue = FaultPlot_3D_Color_SelectedElements(FaultCenter,FaultLengthStrike, FaultLengthDip,
#     FaultStrikeAngle, FaultDipAngle, FaultRakeAngle, PlotInput, 
#     PlotRotation, MinMax_Axis, ColorMinMax, Transparent, SelectedElements)

