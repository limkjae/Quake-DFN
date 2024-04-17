

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


function BuildInputFromBulkGeometry_H()


    InputBulkFileName="Input_BulkFaultGeometry.txt"
    OutputFileName="Input_Discretized_H.jld2"


    ArrangePoint = [0,-3000,1000]
    HowManyDivisionEachLevel = 2
    TotalHierarchyLevel = 8
    MinimumElementsToCut = 30
    ElementPartRoughCount = 2000
    DistDiamRatioCrit = 1

    ############# Plots? ##############
    PlotHMat = 1
    PlotBlock3D = 1
    ###################################


    
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
    
    Switch_StrikeSlip_or_ReverseNormal = Input_Bulk[2,1] 
    ShearModulus = Input_Bulk[2,2]
    PoissonRatio = Input_Bulk[2,3]
    RockDensity = Input_Bulk[2,4]
    DropCrit= Input_Bulk[2,5]
    DropCritNormalStressMultiplier= Input_Bulk[2,6]
    MinimumNS=Input_Bulk[2,7]

    Input_Bulk=Input_Bulk[4:end,:]
    Input_Bulk=Input_Bulk[sortperm(Input_Bulk[:, 17]), :]
    

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
    ##  1.Ctr_X     2.Ctr_Y 3.Ctr_Z 4.St_L	    5.Dip_L	    6.StAng	    7.DipAng	8.LR/RN
    ##  9.a         10.b	11.Dc	12.Theta_i	13. V_i     14. Friction_i 15.NormalStress  
    ##  16. V_Const 17. Bulk Index     18. Bulk Strike Length      19. Bulk Dip Length
    Input_Segment = BulkToSegment(Input_Bulk);
    FaultCount=   size(Input_Segment,1)
    Input_Segment = [Input_Segment ones(FaultCount) zeros(FaultCount)]



    ################################# Group and Sort #################################
    Block_Ctr_Diam, Block_Range_Level, Input_Segment = 
        GroupAndSort_AllLevel(HowManyDivisionEachLevel, TotalHierarchyLevel, MinimumElementsToCut,
                            ArrangePoint, FaultCount, Input_Segment)
    println("Grouping and Sorting Done")




    ############################## Build Hierarchy ####################################
    Ranks, ElementRange_SR = BuildHierarchy(Block_Range_Level, Block_Ctr_Diam, DistDiamRatioCrit, TotalHierarchyLevel) 



    ################################################################################
    #################################### 3D Plot ###################################

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
    ########^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^########
    ################################################################################



    ################################################################################
    ########################### Plot Hierarchical Matrix ###########################
    if PlotHMat == 1
        figure(2)
        println("Plotting Hmatrix Structure")
        clf()
        ax = gca()
        # ax[:set_aspect]("equal")
        for i=1:length(ElementRange_SR[:,1])
            if Ranks[i] > 6 
                c = PyObject(patches.Rectangle((ElementRange_SR[i,3], -ElementRange_SR[i,1]), ElementRange_SR[i,4] - ElementRange_SR[i,3], 
                            -ElementRange_SR[i,2] + ElementRange_SR[i,1], linewidth=1, edgecolor="k", facecolor=[0.2 0.3 1]))                                
            elseif Ranks[i] > 0                 
                c = PyObject(patches.Rectangle((ElementRange_SR[i,3], -ElementRange_SR[i,1]), ElementRange_SR[i,4] - ElementRange_SR[i,3], 
                            -ElementRange_SR[i,2] + ElementRange_SR[i,1], linewidth=1, edgecolor="k", facecolor=[0.5 0.8 1]))
            else           
                c = PyObject(patches.Rectangle((ElementRange_SR[i,3], -ElementRange_SR[i,1]), ElementRange_SR[i,4] - ElementRange_SR[i,3], 
                            -ElementRange_SR[i,2] + ElementRange_SR[i,1], linewidth=1, edgecolor="k", facecolor=[1 0.8 0.4]))
            end
            ax.text( (ElementRange_SR[i,3] + ElementRange_SR[i,4])/2, -(ElementRange_SR[i,1] + ElementRange_SR[i,2])/2, i ,size=8, horizontalalignment="center", verticalalignment="center", color="k") 
            ax.add_patch(c) 
        end
        xlim(0,FaultCount)
        ylim(-FaultCount,0)
    end
    ########^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^########
    ################################################################################






    ########################## Build Stiffness Matrix ##############################

    StiffnessMatrixShearOriginal=zeros(FaultCount,FaultCount)
    StiffnessMatrixNormalOriginal=zeros(FaultCount,FaultCount)

    if Switch_StrikeSlip_or_ReverseNormal == 1
        println("Preparing for discretization")
        StiffnessMatrixShearOriginal, StiffnessMatrixNormalOriginal = 
            BuildMatrixByPartsShear(FaultCount, ElementPartRoughCount, Input_Segment,  ShearModulus, PoissonRatio)
        
    elseif Switch_StrikeSlip_or_ReverseNormal == 2
        println("Preparing for discretization")
        StiffnessMatrixShearOriginal, StiffnessMatrixNormalOriginal = 
            BuildMatrixByPartsNormal(FaultCount, ElementPartRoughCount, Input_Segment,  ShearModulus, PoissonRatio)

    else
        println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        println("!!!!!!! Switch_StrikeSlip_or_ReverseNormal should be either 1 or 2  !!!!!!!")
        println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    end

    ElasticLoadingShearMatrix=zeros(FaultCount,FaultCount)
    for i=1:FaultCount
        for j=1:FaultCount  
            if i==j
                ElasticLoadingShearMatrix[i,j]=0.0;
            else
                ElasticLoadingShearMatrix[i,j]=StiffnessMatrixShearOriginal[i,j]/StiffnessMatrixShearOriginal[i,i];
                # ElasticLoadingNormal[i,j]=StiffnessMatrixNormal[i,j]/StiffnessMatrixShear[i,i];
            end
        end
    end

    ##################### Approximate #####################
    BlockCount = length(ElementRange_SR[:,1])
    ShearStiffness_H = Any[0]
    BlockIndex = 0
    # Ranks = Ranks .* 2
    for i=1:BlockCount
        BlockIndex = BlockIndex + 1
        if Ranks[BlockIndex] > 0
            # push!(ABlock_H,psvdfact(ShearMatrix[ElementRange_SR[i,1]:ElementRange_SR[i,2],ElementRange_SR[i,3]:ElementRange_SR[i,4]], rank=Ranks[i]))
            push!(ShearStiffness_H,pqrfact(ElasticLoadingShearMatrix[ElementRange_SR[i,1]:ElementRange_SR[i,2],ElementRange_SR[i,3]:ElementRange_SR[i,4]], rank=Ranks[i]))
        else 
            push!(ShearStiffness_H,ElasticLoadingShearMatrix[ElementRange_SR[i,1]:ElementRange_SR[i,2],ElementRange_SR[i,3]:ElementRange_SR[i,4]])
        end

    end
    ShearStiffness_H = ShearStiffness_H[2:end]








    println("Stiffness NaN count: ",sum(isnan.(StiffnessMatrixShearOriginal)))

    ########################## Remove Unstable Faults ##############################    
    ReducedStiffnessMatrixShear, ReducedStiffnessMatrixNormal, ReducedInput_Segment=
    CheckTooClose(StiffnessMatrixShearOriginal, StiffnessMatrixNormalOriginal, Input_Segment, Input_Bulk, DropCrit, DropCritNormalStressMultiplier);


    ################################### Save Files #################################
    jldsave("HmatSave.jld2"; ShearStiffness_H, ElasticLoadingShearMatrix, ReducedStiffnessMatrixShear, ReducedStiffnessMatrixNormal, Ranks, ElementRange_SR)

    SaveResults_H(ReducedStiffnessMatrixShear, ReducedStiffnessMatrixNormal, ReducedInput_Segment,
         OutputFileName, ShearModulus, PoissonRatio, RockDensity, Switch_StrikeSlip_or_ReverseNormal, MinimumNS,
         Ranks, ElementRange_SR, ShearStiffness_H);

end




Input = BuildInputFromBulkGeometry_H()