
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

# figure(10); clf(); PyPlot.plot(log10.(ResultV[:,1:30:end])); xlabel("Record Step")
PlotStep = 135

PlotRotation = [35,-30]
Transparent = 0 # 1 for transparent fault plot. 0 for no-transparency
Edge = 0 # 0 for no element boudary. 1 for plotting element boundary
MinMax_Axis = 0 # 0 for automatically selected axis minimim and maximum 
# MinMax_Axis=[-2000 2000; -2000 2000; -4000 0]
LoadingFaultPlot = 0 # 1 to plot constant velocity faults. 

ShowDay = 0 # If 1, day is shown in the location 
DayLocation = [0,0,1000]


################## What to Plot ? ###################
PlotInput=log10.(ResultV[PlotStep,:]); ColorMinMax=[-12,0] 
# PlotInput= Result_NormalStress[PlotStep,:] ; ColorMinMax=0
# PlotInput=log10.(ResultPressure[PlotStep,:]); ColorMinMax=[3,6]
# PlotInput= Result_NormalStress[PlotStep,:] -  Fault_NormalStress; ColorMinMax=[-1e6,1e6]
# PlotInput=ResultDisp[PlotStep,:]; ColorMinMax=0 
#############---------------------------#############


############# Saving Multiple Figures ###############
Animation_Save = 0 # 1 for save
StepBegin = 10 # first record step
StepEnd = 100
StepInterval = 10
###################^^^^^^^^^^^^^^^###################

##############+++++++++++++++++++++++++++++++++++++++++++++++++++++++++################
#######################################################################################



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
Fault_BulkIndex= load(FileNameInput, "Fault_BulkIndex")
FaultMass= load(FileNameInput, "FaultMass")

if LoadingFaultPlot==1
LoadingFaultCount=0
end

figure(1)
clf()
MaxVaule, MinValue = FaultPlot_3D_Color_General(FaultCenter,FaultLengthStrike, FaultLengthDip,
    FaultStrikeAngle, FaultDipAngle, FaultLLRR, PlotInput, 
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




if Animation_Save == 1

    # PlotStepI=PlotStep
    isdir("3DPlot") || mkdir("3DPlot")
    for i=StepBegin:StepInterval:StepEnd
        PlotStep = i#FaultPlot_3D(FaultCenter[1:FaultCount-2,:],FaultLengthStrike[1:FaultCount-2], FaultLengthDip[1:FaultCount-2], FaultStrikeAngle[1:FaultCount-2], FaultDipAngle[1:FaultCount-2], FaultLLRR[1:FaultCount-2])
        #FaultPlot_3D_ColorDisp(FaultCenter[1:FaultCount-2,:],FaultLengthStrike[1:FaultCount-2], FaultLengthDip[1:FaultCount-2], FaultStrikeAngle[1:FaultCount-2], FaultDipAngle[1:FaultCount-2], FaultLLRR[1:FaultCount-2],ResultV[:,1:FaultCount-2], PlotStep)
        clf()
        PlotInput=log10.(ResultV[PlotStep,:])

        MaxVaule, MinValue = FaultPlot_3D_Color_General(FaultCenter,FaultLengthStrike, FaultLengthDip,
        FaultStrikeAngle, FaultDipAngle, FaultLLRR, PlotInput, 
        PlotRotation, MinMax_Axis, ColorMinMax, Transparent, Edge, LoadingFaultCount)
        figure(1).canvas.draw()
        ax = subplot(projection="3d")
        PlotTime=ResultTime[PlotStep]/60/60/24
        if ShowDay ==1 
        # ax.text(-5000, 10, 1300, "Day: ",size=10)
        ax.text(DayLocation[1], DayLocation[2], DayLocation[3], PlotTime,size=10)    
        end    
        figure(1).canvas.draw()

        println(PlotStep)
        
        PyPlot.savefig("3DPlot/" * ResultName * "_" * string(i) * ".png")
    end

end



# #Single element locator
# # figure(1)
# SelectedElements = [100]
# # clf()
# MaxVaule, MinValue = FaultPlot_3D_Color_SelectedElements(FaultCenter,FaultLengthStrike, FaultLengthDip,
#     FaultStrikeAngle, FaultDipAngle, FaultLLRR, PlotInput, 
#     PlotRotation, MinMax_Axis, ColorMinMax, Transparent, SelectedElements)

 