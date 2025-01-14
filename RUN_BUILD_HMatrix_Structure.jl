


################################################################################
###### This file generate Hmatrix Structure File from Bulk Fault Geometry ######
###### H-matrix approximation should  be 

using DelimitedFiles
using Base
using PyPlot
using PyCall
using JLD2
using LowRankApprox
using Clustering
using LinearAlgebra
@pyimport matplotlib.patches as patches

pygui(true)

include("Functions_BuildInputFile.jl")
include("Functions_OKADA3D.jl")
include("Results/Functions_Plot.jl")
include("Functions_Hmatrix.jl")

function BuildHMatStructure()

    InputBulkFileName="Input_BulkFaultGeometry.txt"
    OutputFileName="Input_HmatrixStructure.jld2"

    ##########################################################################
    #####----- Hmatrix compression detail ----#####
    TotalHierarchyLevel = 8
    MinimumElementsToCut = 30
    DistDiamRatioCrit = 1.0     
    ArrangePoint = [10000,6000,0] ### Arrange point
    # ArrangePoint = [0,0,-1000] ### Arrange point

    #####---------   HMatrix Plots?  --------#####
    PlotHMat = 1 # HMatrix structure plot
    PlotBlock3D = 0
    ##########################################################################
    

    ElementPartRoughCount = 2000
    
    ############################# Read Bulk Input ##################################
    ######++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++######
    ## Bulk File Order 
    ##  1.SwitchSSRN  2.ShearMod  3.PoiRat  4.R_Density   
    ##  5. Crit_TooClose   6. TooCloseNormal_Multiplier
    ##  ----------------------------------------------------------------------------
    ##  1.Ctr_X     2.Ctr_Y 3.Ctr_Z 4.St_L	    5.Dip_L	    6.StAng	    7.DipAng	8.LR/RN
    ##  9.a         10.b	11.Dc	12.Theta_i	13. V_i     14. Friction_i 15.NormalStress at surface [Pa]  
    ##  16. NoarmalStress Gradient [Pa] 17. V_Const     18. Minimum Segment Length

    Input_Bulk=readdlm(InputBulkFileName)
    Input_BulkOriginal=copy(Input_Bulk)

    
    Switch_StrikeSlip_or_ReverseNormal = Input_Bulk[2,1] 
    ShearModulus = Input_Bulk[2,2]
    PoissonRatio = Input_Bulk[2,3]
    RockDensity = Input_Bulk[2,4]
    DropCrit= Input_Bulk[2,5]
    DropCritNormalStressMultiplier= Input_Bulk[2,6]
    MinimumNS=Input_Bulk[2,7]

    Input_Bulk=Input_Bulk[4:end,:]
    Input_Bulk=Input_Bulk[sortperm(Input_Bulk[:, 17]), :]

    BulkFaultcount = length(Input_Bulk[:,1])
    # Adjust LRRN to rake angle
    if Switch_StrikeSlip_or_ReverseNormal == 1
        for BulkIndex = 1: BulkFaultcount
            if Input_Bulk[BulkIndex,8] == -1.0
                Input_Bulk[BulkIndex,8] = 0.0
            else        
                Input_Bulk[BulkIndex,8] = 180.0
            end    
        end
    elseif Switch_StrikeSlip_or_ReverseNormal ==2
        for BulkIndex = 1: BulkFaultcount
            if Input_Bulk[BulkIndex,7] < 90.0
                if Input_Bulk[BulkIndex,8] == -1.0
                    Input_Bulk[BulkIndex,8] = 90.0
                else        
                    Input_Bulk[BulkIndex,8] = 270.0
                end    
            else 
                if Input_Bulk[BulkIndex,8] == -1.0
                    Input_Bulk[BulkIndex,8] = 270.0
                else        
                    Input_Bulk[BulkIndex,8] = 90.0
                end    
            end
        end
    end   


    # Adjust if positive depth exists
    for i in eachindex(Input_Bulk[:,1])
        if Input_Bulk[i,3]<Input_Bulk[i,5]/2*sind(Input_Bulk[i,7])
            println("Caution! Fault ",i," may have negative depth")
        end
    end

    ########^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^########
    ################################################################################



    ############################### Bulk to Segment ################################
    ######++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++######
    ## Segment File Order 
    ##  1.Ctr_X     2.Ctr_Y 3.Ctr_Z 4.St_L	    5.Dip_L	    6.StAng	    7.DipAng	8.RakeAngle
    ##  9.a         10.b	11.Dc	12.Theta_i	13. V_i     14. Friction_i 15.NormalStress  
    ##  16. V_Const 17. Bulk Index     18. Bulk Strike Length      19. Bulk Dip Length
    Input_Segment = BulkToSegment(Input_Bulk);
    FaultCount=   size(Input_Segment,1)
    Input_Segment = [Input_Segment ones(FaultCount) zeros(FaultCount)]
    Input_Segment = Input_Segment[sortperm(Input_Segment[:, 16]), :] # move the loading faults to the top
    


    #################################  Hmatrix  ######################################


    #####------------------------- Group and Sort -------------------------------#####
    HowManyDivisionEachLevel = 2
    Block_Ctr_Diam, Block_Range_Level, Input_Segment, LoadingFaultExist = 
        GroupAndSort_AllLevel(HowManyDivisionEachLevel, TotalHierarchyLevel, MinimumElementsToCut,
                            ArrangePoint, FaultCount, Input_Segment)
    println("Grouping and Sorting Done")


    #####------------------------- Build Hierarchy -------------------------------#####
    Admissible, ElementRange_SR = BuildHierarchy(Block_Range_Level, Block_Ctr_Diam, DistDiamRatioCrit, TotalHierarchyLevel, LoadingFaultExist) 


    println("total ", length(ElementRange_SR[:,1]), " blocks generated (Low rank:", sum(Admissible)," blocks)")

    #####---------------------------     Plots     -------------------------------#####
    #####                               3D Plot                                   #####
    if PlotBlock3D == 1
        println("Plotting 3D Group")
        PlotRotation=[45,-30]
        Edge = 0
        Transparent = 0
        MinMax_Axis=0 # automatically detect max and min 
        LoadingFaultCount=0 
        ColorMinMax = 0  
        figure(1)
        clf()
        MaxVaule, MinValue = FaultPlot_3D_Color_General(Input_Segment[:,1:3],Input_Segment[:,4], Input_Segment[:,5],
            Input_Segment[:,6], Input_Segment[:,7], Input_Segment[:,8], Input_Segment[:,20], 
            PlotRotation, MinMax_Axis, ColorMinMax, Transparent, Edge, LoadingFaultCount)

        plotforcbar=  scatter([1,1],[1,1],0.1, [MinValue,MaxVaule], cmap="jet")
        colorbar(plotforcbar, pad=0.15)
        figure(1).canvas.draw()
        xlabel("x")
        ylabel("y")
    end

    #####                     Hierarchical Matrix Structure Plot                  #####
    if PlotHMat == 1
        figure(2)
        println("Plotting Hmatrix Structure")
        clf()
        ax = gca()
        # ax[:set_aspect]("equal")
        for i=1:length(ElementRange_SR[:,1])
            if Admissible[i] > 0   
                c = PyObject(patches.Rectangle((ElementRange_SR[i,1]-1, -ElementRange_SR[i,3]+1), ElementRange_SR[i,2] - ElementRange_SR[i,1]+1, 
                            -ElementRange_SR[i,4] + ElementRange_SR[i,3]-1, linewidth=1, edgecolor="k", facecolor=[0.4  0.4  1]))    
            else           
                c = PyObject(patches.Rectangle((ElementRange_SR[i,1]-1, -ElementRange_SR[i,3]+1), ElementRange_SR[i,2] - ElementRange_SR[i,1]+1, 
                            -ElementRange_SR[i,4] + ElementRange_SR[i,3]-1, linewidth=1, edgecolor="k",  facecolor=[1 0.4 0.4]))
            end
            # ax.text( (ElementRange_SR[i,3] + ElementRange_SR[i,4])/2, -(ElementRange_SR[i,1] + ElementRange_SR[i,2])/2, i ,size=8, horizontalalignment="center", verticalalignment="center", color="k") 
            ax.add_patch(c) 
        end
        xlim(0,FaultCount)
        ylim(-FaultCount,0)
    end
    ########^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^########
        
    #################### end of Hmatrix Sorting and Plotting #######################



    ###############################################################################
    ############################### Save Input Files ##############################
    save(OutputFileName, 
    "Block_Ctr_Diam", Block_Ctr_Diam, 
    "Block_Range_Level", Block_Range_Level, 
    "Input_Segment", Input_Segment, 
    "LoadingFaultExist", LoadingFaultExist,
    "HowManyDivisionEachLevel", HowManyDivisionEachLevel, 
    "MinimumElementsToCut", MinimumElementsToCut, 
    "ArrangePoint", ArrangePoint, 
    "Admissible", Admissible, 
    "ElementRange_SR", ElementRange_SR, 
    "Switch_StrikeSlip_or_ReverseNormal", Switch_StrikeSlip_or_ReverseNormal, 
    "ShearModulus", ShearModulus, 
    "PoissonRatio", PoissonRatio, 
    "RockDensity", RockDensity, 
    "DropCrit", DropCrit, 
    "DropCritNormalStressMultiplier", DropCritNormalStressMultiplier, 
    "MinimumNS", MinimumNS,
    "Input_Bulk", Input_BulkOriginal)

    

end


Input = BuildHMatStructure()
