
using DelimitedFiles
using PyPlot
using PyCall

pygui(true)
include("scripts/Functions_BuildInputFile.jl")
include("scripts/Functions_Plot.jl")
InputBulkFileName="Input_BulkFaultGeometry.txt"


function PlotBulk()

    Input_Bulk=readdlm(InputBulkFileName)
    RakeLRRN = Int(Input_Bulk[2,1])
    Input_Bulk=Input_Bulk[4:end,:] 
 
    if size(Input_Bulk, 2) == 18
        println("Rectangle")
        Input_Bulk = Input_Bulk[sortperm(Input_Bulk[:, 17], rev=false), :] # move 
        RorT = "R"
    elseif  size(Input_Bulk, 2) == 20
        println("Triangle")
        Input_Bulk = Input_Bulk[sortperm(Input_Bulk[:, 19], rev=false), :] # move 
        RorT = "T"
    else 
        error("Input Bulk Fault Geometry file should have 18 or 20 columns")
    end
    #######################################################################################
    ################################  Input Rectangle #####################################
    ##  1.Ctr_X     2.Ctr_Y 3.Ctr_Z 4.St_L	    5.Dip_L	    6.StAng	    7.DipAng	8.LR/RN
    ##  9.a         10.b	11.Dc	12.Theta_i	13. V_i     14. Friction_i 15.NormalStress at surface [Pa]  
    ##  16. NoarmalStress Gradient [Pa] 17. V_Const     18. Minimum Segment Length

    ################################  Input Triangle  #####################################
    ####  123. P1xyz    456. P2xyz  789. P3xyz      10.Rake
    ####  11.a          12.b	    13.Dc	        14.Theta_i	15. V_i     16. Friction_i 
    ####  17.NormalStress at surface [Pa]  18. NoarmalStress Gradient [Pa] 
    ####  19. V_Const     20. Minimum Segment Length

    PlotInput = Input_Bulk[:,10]; ColorMinMax=0    
    ##### PlotInput = log10.(Input_Bulk[:,15]); ColorMinMax=[-20,0]    

    ############################### Figure Configuration ##################################
    PlotRotation=[30,-50]
    Transparent = 1 # 1 for transparent fault plot
    Edge = 1 # 0 for no element boudary 
    MinMax_Axis=0
    LoadingFaultPlot = 0 # 1 to plot constant velocity faults. 
    #######################################################################################

    if LoadingFaultPlot == 0
        if RorT == "R"   
            LoadingFaultCount = sum(x->x>0, Input_Bulk[:,17]) 
        else
            LoadingFaultCount = sum(x->x>0, Input_Bulk[:,19]) 
        end
        Input_Bulk = Input_Bulk[1:end-LoadingFaultCount,:]
        PlotInput = PlotInput[1:end-LoadingFaultCount]
    end
    FaultCount = length(Input_Bulk[:,1]) 

    figure(1)
    clf()

    if RorT == "R"     
        MaxVaule, MinValue = FaultPlot_3D_Color_General(Input_Bulk[:,1:3],
            Input_Bulk[:,4], Input_Bulk[:,5], Input_Bulk[:,6], Input_Bulk[:,7], Input_Bulk[:,8], PlotInput, 
            PlotRotation, MinMax_Axis, ColorMinMax, Transparent, Edge, 0)

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

            ax.set_aspect("equal")

    elseif RorT == "T"

        if ColorMinMax == 0 
        MaxValue=maximum(PlotInput)
        MinValue=minimum(PlotInput)
        else
        MaxValue=ColorMinMax[2]
        MinValue=ColorMinMax[1]
        end

        P1, P2, P3 = Input_Bulk[:,1:3], Input_Bulk[:,4:6], Input_Bulk[:,7:9]
        figure(1)
        fig = figure(1)
        clf()
        art3d = PyObject(PyPlot.art3D)
        ax = subplot(projection="3d")
        if Edge == 0 
            edge_color = [0.2, 0.2, 0.2, 0.0]
        else
            edge_color = [0.2, 0.2, 0.2, 0.2]
        end

        for ElemIdx = 1:FaultCount
            cm = get_cmap(:jet)
            PlotValue=(PlotInput[ElemIdx]-MinValue)/(MaxValue-MinValue)

            if Transparent ==0
                face_color = [cm(PlotValue)[1], cm(PlotValue)[2],cm(PlotValue)[3],1.0]
            else
                face_color = [cm(PlotValue)[1], cm(PlotValue)[2],cm(PlotValue)[3],0.5]
            end

            verts = ((P1[ElemIdx,:],P2[ElemIdx,:],P3[ElemIdx,:]), )
            p3c = PyObject(art3d.Poly3DCollection(verts))
            pycall(ax.add_collection3d, PyAny, p3c)
            pycall(p3c.set_facecolor, PyAny, face_color)
            pycall(p3c.set_edgecolor, PyAny, edge_color)
            ax.view_init(PlotRotation[1], PlotRotation[2])
        end

        plotforcbar=  scatter([P1[1,1],P1[2,1]],[P1[1,2],P1[2,2]],0.1, [MinValue,MaxValue], cmap="jet")
        colorbar(plotforcbar, pad=0.15)
        figure(1).canvas.draw()
        ax.set_aspect("equal")
        xlabel("x")
        ylabel("y")
        zlabel("z")
    end



end

PlotBulk()

