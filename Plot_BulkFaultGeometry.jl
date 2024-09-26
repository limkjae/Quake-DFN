
using DelimitedFiles
using PyPlot
using PyCall

pygui(true)
include("Functions_BuildInputFile.jl")
include("Results/Functions_Plot.jl")
InputBulkFileName="Input_BulkFaultGeometry.txt"
Input_Bulk=readdlm(InputBulkFileName)
Input_Bulk=Input_Bulk[4:end,:]


#######################################################################################
############################### Figure Configuration ##################################

PlotRotation=[30,-50]
Transparent = 1 # 1 for transparent fault plot
Edge = 1 # 0 for no element boudary 
MinMax_Axis=0
LoadingFaultCount=0 # 1 to plot constant velocity faults. 
Input_Bulk = Input_Bulk[1:end,:]
PlotInput=Input_Bulk[:,8]; ColorMinMax=0 
# PlotInput=log10.(Input_Bulk[:,13]); ColorMinMax=0

##### Input Code
##### 1.Ctr_X     2.Ctr_Y 3.Ctr_Z 4.St_L	    5.Dip_L	    6.StAng	    7.DipAng	8.LR/RN
#####  9.a         10.b	    11.Dc	12.Theta_i	13. V_i     14. Friction_i 15.NormalStress at surface [Pa]  
#####  16. NoarmalStress Gradient [Pa] 17. V_Const     18. Minimum Segment Length


##############+++++++++++++++++++++++++++++++++++++++++++++++++++++++++################
#######################################################################################


figure(1)
clf()
MaxVaule, MinValue = FaultPlot_3D_Color_General(Input_Bulk[:,1:3],
    Input_Bulk[:,4], Input_Bulk[:,5], Input_Bulk[:,6], Input_Bulk[:,7], Input_Bulk[:,8], PlotInput, 
    PlotRotation, MinMax_Axis, ColorMinMax, Transparent, Edge, LoadingFaultCount)
    ax = subplot(projection="3d")
    xlabel("x")
    ylabel("y")
plotforcbar=  scatter([1,1],[1,1],0.1, [MinValue,MaxVaule], cmap="jet")
colorbar(plotforcbar, pad=0.15)
figure(1).canvas.draw()
