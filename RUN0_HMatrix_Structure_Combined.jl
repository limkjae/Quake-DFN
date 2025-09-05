
################################################################################
###### This file generate Hmatrix Structure File from Bulk Fault Geometry ######

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
    MinimumElementsToCut = 10
    DistDiamRatioCrit = 1.0     
    ArrangePoint = [10000,6000,10000] ### Arrange point
    # ArrangePoint = [0,0,-1000] ### Arrange point

    #####---------   HMatrix Plots?  --------#####
    PlotHMat = 1 # HMatrix structure plot
    PlotBlock3D = 0
    ##########################################################################


    ######################## Read Bulk Input File ########################
    Input_Segment, LoadingFaultCount, ShearModulus, PoissonRatio, RockDensity, 
    Switch_StrikeSlip_or_ReverseNormal, DropCrit, DropCritNormalStressMultiplier, MinimumNS, RorT, FaultCount =
        ReadBulkInput(InputBulkFileName)


    Input_Segment = [Input_Segment ones(FaultCount) zeros(FaultCount)]


    ######################## Build Hmatrix structure ########################
    Block_Ctr_Diam, Block_Range_Level, Input_Segment, LoadingFaultExist = 
        GroupAndSort_AllLevel(TotalHierarchyLevel, MinimumElementsToCut,
                            ArrangePoint, Input_Segment, RorT)  
    println("Grouping and Sorting Done")

    Admissible, ElementRange_SR = BuildHierarchy(Block_Range_Level, Block_Ctr_Diam, 
                                DistDiamRatioCrit, TotalHierarchyLevel, LoadingFaultExist) 

    println("total ", length(ElementRange_SR[:,1]), " blocks generated (Low rank:", sum(Admissible)," blocks)")




    ###############################################################################
    ############################### Save Input Files ##############################
    save(OutputFileName, 
    "LoadingFaultCount", LoadingFaultCount,
    "Block_Ctr_Diam", Block_Ctr_Diam, 
    "Block_Range_Level", Block_Range_Level, 
    "Input_Segment", Input_Segment, 
    "LoadingFaultExist", LoadingFaultExist,
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
    "RorT", RorT)
    println("Input files saved to ", OutputFileName) 

     



  ############### Plots ####################


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




    # if RorT == "R"
    #     PlotInput = Input_Segment[:,21]; ColorMinMax= 0
    #     MaxValue=maximum(PlotInput)
    #     MinValue=minimum(PlotInput)
        
    #     figure(1)
    #     fig = figure(1)
    #     clf()
    #     art3d = PyObject(PyPlot.art3D)
    #     ax = subplot(projection="3d")
    #     for ElemIdx = 1:length(PlotInput) - LoadingFaultCount
            
    #         cm = get_cmap(:prism)
    #         PlotValue=(PlotInput[ElemIdx]-MinValue)/(MaxValue-MinValue)
    #         FaultStrikeAngle = Input_Segment[:,6]
    #         FaultDipAngle = Input_Segment[:,7]
    #         FaultLengthStrike = Input_Segment[:,4]
    #         FaultLengthDip = Input_Segment[:,5]
    #         FaultCenter = Input_Segment[:,1:3]
    #         face_color = [cm(PlotValue)[1], cm(PlotValue)[2],cm(PlotValue)[3],1.0]


    #         RotMatStrike=[cosd(FaultStrikeAngle[ElemIdx]) -sind(FaultStrikeAngle[ElemIdx]) 0
    #             sind(FaultStrikeAngle[ElemIdx]) cosd(FaultStrikeAngle[ElemIdx]) 0
    #             0  0 1]
    #         RotMatDip=[1 0  0
    #         0 cosd(FaultDipAngle[ElemIdx]) -sind(FaultDipAngle[ElemIdx])
    #         0 sind(FaultDipAngle[ElemIdx]) cosd(FaultDipAngle[ElemIdx])]
                    
    #         p1=RotMatStrike*RotMatDip*[FaultLengthStrike[ElemIdx]/2;-FaultLengthDip[ElemIdx]/2;0] + [FaultCenter[ElemIdx,1]; FaultCenter[ElemIdx,2]; -FaultCenter[ElemIdx,3]];
    #         p2=RotMatStrike*RotMatDip*[-FaultLengthStrike[ElemIdx]/2;-FaultLengthDip[ElemIdx]/2;0] + [FaultCenter[ElemIdx,1]; FaultCenter[ElemIdx,2]; -FaultCenter[ElemIdx,3]];
    #         p3=RotMatStrike*RotMatDip*[-FaultLengthStrike[ElemIdx]/2;+FaultLengthDip[ElemIdx]/2;0]+ [FaultCenter[ElemIdx,1]; FaultCenter[ElemIdx,2]; -FaultCenter[ElemIdx,3]];
    #         p4=RotMatStrike*RotMatDip*[FaultLengthStrike[ElemIdx]/2;+FaultLengthDip[ElemIdx]/2;0]+ [FaultCenter[ElemIdx,1]; FaultCenter[ElemIdx,2]; -FaultCenter[ElemIdx,3]];
            
    #         # verts2 = ([tuple(p1...); tuple(p2...); tuple(p3...); tuple(p4...)],)

    #         verts =([tuple(p1...); tuple(p2...); tuple(p3...); tuple(p4...)],)
    #         p3c = PyObject(art3d.Poly3DCollection(verts))
    #         pycall(ax.add_collection3d, PyAny, p3c)

    #         # face_color = [0.3, 0.8, 0.3, 0.5]         
    #         edge_color = [0.2, 0.2, 0.2, 0.0]

    #         pycall(p3c.set_facecolor, PyAny, face_color)
    #         pycall(p3c.set_edgecolor, PyAny, edge_color)
    #         ax.view_init(17, -70)
    #     end
    # else
    #     PlotInput = Input_Segment[:,21]; ColorMinMax= 0
    #     MaxValue=maximum(PlotInput)
    #     MinValue=minimum(PlotInput)
        
    #     figure(1)
    #     fig = figure(1)
    #     clf()
    #     art3d = PyObject(PyPlot.art3D)
    #     ax = subplot(projection="3d")
    #     for ElemIdx = 1:FaultCount- LoadingFaultCount
    #         cm = get_cmap(:prism)
    #         P1 = Input_Segment[:,1:3]
    #         P2 = Input_Segment[:,4:6]
    #         P3 = Input_Segment[:,7:9]
    #         PlotValue=(PlotInput[ElemIdx]-MinValue)/(MaxValue-MinValue)

    #         face_color = [cm(PlotValue)[1], cm(PlotValue)[2],cm(PlotValue)[3],1.0]

    #         verts = ((P1[ElemIdx,:],P2[ElemIdx,:],P3[ElemIdx,:]), )
    #         p3c = PyObject(art3d.Poly3DCollection(verts))
    #         pycall(ax.add_collection3d, PyAny, p3c)

    #         # face_color = [0.3, 0.8, 0.3, 0.5]         
    #         edge_color = [0.2, 0.2, 0.2, 0.0]

    #         pycall(p3c.set_facecolor, PyAny, face_color)
    #         pycall(p3c.set_edgecolor, PyAny, edge_color)
    #         ax.view_init(17, -70)
    #     end
    # end

end

BuildHMatStructure()
