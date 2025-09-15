
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




end

BuildHMatStructure()
