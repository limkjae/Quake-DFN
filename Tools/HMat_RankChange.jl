
using PyPlot
using PyCall
using DelimitedFiles
using JLD2
using LinearAlgebra
using Printf
using SpecialFunctions
using StaticArrays
using LowRankApprox
using Distributed
using Statistics
@pyimport matplotlib.patches as patches
pygui(true)

#########################################################################
####### This script compress the HMatrix with different Tolerance #######
#######               Original Matrix Should exists               #######
#######         Rank can be changed by adjusting Tolerance        #######
#######             HMatrix Structure Cannot be changed           #######
#########################################################################



######################## Recompression Tolerance ##########################
Tolerance = 1e3 # pascal for 1m slip (More approximaion for higher Tolerance)
PlotHMat = 1 # HMatrix structure plot
###########################################################################

include("../Functions_Solvers.jl")
include("../Functions_RSFDFN3DMain_H.jl")
include("../Results/Functions_Plot.jl")
include("../QuickParameterAdjust.jl")
include("../Functions_Hmatrix.jl")
LoadingInputFileName="Input_Discretized.jld2" 


function rank_change()

    # jldsave("HmatSave.jld2"; StiffnessMatrixShearOriginal, StiffnessMatrixNormalOriginal, Ranks, ElementRange_SR)
    Admissible = load(LoadingInputFileName,"Admissible")
    FaultCount= load(LoadingInputFileName, "FaultCount")
    Ranks_Shear= load(LoadingInputFileName, "Ranks_Shear") # figure(11); plot(Ranks)
    Ranks_Normal= load(LoadingInputFileName, "Ranks_Normal") # figure(11); plot(Ranks)
    ElementRange_SR = load(LoadingInputFileName, "ElementRange_SR")
    ShearStiffness_H = load(LoadingInputFileName, "ShearStiffness_H")
    NormalStiffness_H = load(LoadingInputFileName, "NormalStiffness_H")
    NormalStiffnessZero = load(LoadingInputFileName, "NormalStiffnessZero")
    StiffnessMatrixShear= load(LoadingInputFileName, "StiffnessMatrixShear")
    StiffnessMatrixNormal= load(LoadingInputFileName, "StiffnessMatrixNormal")




    figure(11)
    plot(Ranks_Shear,"b")
    figure(12)
    plot(Ranks_Normal,"b")



    ###################################################################################
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
    #########^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^##########
    ###################################################################################
        

    ###################################################################################
    ############################### HMatrix Approximate ###############################
        BlockCount = length(ElementRange_SR[:,1])
        ShearStiffness_H = Any[0]
        Ranks_Shear = zeros(Int, BlockCount)
        NormalStiffness_H = Any[0]
        Ranks_Normal = zeros(Int, BlockCount)
        BlockIndex = 0
        println("compressing")
        for i=1:BlockCount
            BlockIndex = BlockIndex + 1
            
            if Admissible[BlockIndex] > 0
                
                OrigianlMatrixToApproximate = StiffnessMatrixShear[ElementRange_SR[i,3]:ElementRange_SR[i,4],ElementRange_SR[i,1]:ElementRange_SR[i,2]]
                ApproxMatrixS = pqrfact(OrigianlMatrixToApproximate, atol = Tolerance)
                push!(ShearStiffness_H,ApproxMatrixS)
                Ranks_Shear[BlockIndex] = size(ApproxMatrixS[:Q],2)
                
                OrigianlMatrixToApproximate = StiffnessMatrixNormal[ElementRange_SR[i,3]:ElementRange_SR[i,4],ElementRange_SR[i,1]:ElementRange_SR[i,2]]
                ApproxMatrixN = pqrfact(OrigianlMatrixToApproximate, atol = Tolerance)
                push!(NormalStiffness_H,ApproxMatrixN)
                Ranks_Normal[BlockIndex] = size(ApproxMatrixN[:Q],2)
            else 
                push!(ShearStiffness_H,StiffnessMatrixShear[ElementRange_SR[i,3]:ElementRange_SR[i,4],ElementRange_SR[i,1]:ElementRange_SR[i,2]])
                push!(NormalStiffness_H,StiffnessMatrixNormal[ElementRange_SR[i,3]:ElementRange_SR[i,4],ElementRange_SR[i,1]:ElementRange_SR[i,2]])
            end

        end
        ShearStiffness_H = ShearStiffness_H[2:end]
        NormalStiffness_H = NormalStiffness_H[2:end]

    ########################## end of HMatrix Approximation ##########################

    figure(11)
    title("Ranks in Shear H-Matrix (blue: original, Red: changed)")
    plot(Ranks_Shear,"r")
    figure(12)
    title("Ranks in Normal H-Matrix (blue: original, Red: changed)")
    plot(Ranks_Normal,"r")


    
    file = jldopen(LoadingInputFileName, "a+")
    Base.delete!(file, "Ranks_Shear") 
    Base.delete!(file, "Ranks_Normal") 
    Base.delete!(file, "ShearStiffness_H") 
    Base.delete!(file, "NormalStiffness_H") 


    write(file, "Ranks_Shear", Ranks_Shear) 
    write(file, "Ranks_Normal", Ranks_Normal) 
    write(file, "ShearStiffness_H", ShearStiffness_H) 
    write(file, "NormalStiffness_H", NormalStiffness_H) 


    close(file)

end

rank_change()