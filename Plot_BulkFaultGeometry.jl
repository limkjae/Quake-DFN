
using DelimitedFiles
using PyPlot
using PyCall

pygui(true)
include("Functions_BuildInputFile.jl")
include("Results/Functions_Plot.jl")
InputBulkFileName="Input_BulkFaultGeometry.txt"


function PlotBulk()

    Input_Bulk=readdlm(InputBulkFileName)
    RakeLRRN = Int(Input_Bulk[2,1])
    Input_Bulk=Input_Bulk[4:end,:]
    Input_Bulk = Input_Bulk[sortperm(Input_Bulk[:, 17], rev=false), :] # move the loading faults to the end

    #######################################################################################
    ####################################  Input  ##########################################
    ##  1.Ctr_X     2.Ctr_Y 3.Ctr_Z 4.St_L	    5.Dip_L	    6.StAng	    7.DipAng	8.LR/RN
    ##  9.a         10.b	11.Dc	12.Theta_i	13. V_i     14. Friction_i 15.NormalStress at surface [Pa]  
    ##  16. NoarmalStress Gradient [Pa] 17. V_Const     18. Minimum Segment Length
    PlotInput = Input_Bulk[:,13]; ColorMinMax=0    
    # PlotInput = log10.(Input_Bulk[:,13]); ColorMinMax=[-30,0]    

    ############################### Figure Configuration ##################################
    PlotRotation=[30,-50]
    Transparent = 1 # 1 for transparent fault plot
    Edge = 1 # 0 for no element boudary 
    MinMax_Axis=0
    LoadingFaultPlot = 1 # 1 to plot constant velocity faults. 
    #######################################################################################

    if LoadingFaultPlot == 0
        LoadingFaultCount = sum(x->x>0, Input_Bulk[:,17]) 
        Input_Bulk = Input_Bulk[1:end-LoadingFaultCount,:]
    end



    figure(1)
    clf()
    MaxVaule, MinValue = FaultPlot_3D_Color_General(Input_Bulk[:,1:3],
        Input_Bulk[:,4], Input_Bulk[:,5], Input_Bulk[:,6], Input_Bulk[:,7], Input_Bulk[:,8], PlotInput, 
        PlotRotation, MinMax_Axis, ColorMinMax, Transparent, Edge, 0)
        ax = subplot(projection="3d")
        xlabel("x")
        ylabel("y")
    plotforcbar=  scatter([1,1],[1,1],0.1, [MinValue,MaxVaule], cmap="jet")
    cbar  = colorbar(plotforcbar, pad=0.15)
    # cbar.set_ticks([0, 0.25, 0.5, 0.75,1])
    # cbar.set_ticklabels(["L", "Reverse", "Right Lateral", "Normal", "Left Lateral"])
    figure(1).canvas.draw()


    PlotBulk_SenseOfSlip(RakeLRRN, Input_Bulk, PlotRotation, Transparent, Edge, MinMax_Axis)
 
end

PlotBulk()