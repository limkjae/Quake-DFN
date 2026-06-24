using DelimitedFiles
using Base
# using PyPlot
# using PyCall
using GLMakie
using Statistics
using JLD2
using LowRankApprox
using Clustering
using LinearAlgebra
using Base.Threads
# @pyimport matplotlib.patches as patches
# pygui(true)


include("../scripts/Functions_H2.jl") 
# ElementCount = 4800

H2File = "Input_Discretized.jld2"


FaultCenter, Block_Range_Level, Admissible, ElementRange_SR, Rank, BlockToClusterSR,
FaultCount, Block_Ctr_Diam, BlockCount, ClusterCount, IsLeaf, IsSingleMom, 
DaughterElem, Cluster_ClusterRange, NonLeafClusterIdx, LeafClusterIdx, SingleMomClusterIdx, RegularClusterIdx, 
Cluster_SkeletonGlobal_S, K_d_S, Cluster_W_S,
Cluster_SkeletonGlobal_N, K_d_N, Cluster_W_N =
load(H2File, "FaultCenter", "Block_Range_Level","Admissible","ElementRange_SR", "Rank","BlockToClusterSR",
    "FaultCount", "Block_Ctr_Diam", "BlockCount", "ClusterCount", "IsLeaf", "IsSingleMom",
    "DaughterElem", "Cluster_ClusterRange", "NonLeafClusterIdx",  "LeafClusterIdx", "SingleMomClusterIdx", "RegularClusterIdx",
    "Cluster_SkeletonGlobal_S","K_d_S",  "Cluster_W_S",  
    "Cluster_SkeletonGlobal_N", "K_d_N", "Cluster_W_N")

println("H2Matrix Loaded")



######################## Let's solve #######################

ThreadCount = 5

Disp = ones(FaultCount) 
# Disp = rand(length(StiffnessMatrixShear[1,:])) 


StressShear = zeros(FaultCount)
StressNormal = zeros(FaultCount)



AdmissibleBlockSep, InAdmissibleBlockSep = 
            Ad_InAd_Blocks_EachThreads(BlockCount, Admissible, ThreadCount, FaultCount, 
                                Block_Range_Level, BlockToClusterSR)


    Stress = zeros(FaultCount)
    D_Cluster = [zeros(Float64, Rank) for _ in 1:ClusterCount]
    Sig_Cluster = [zeros(Float64, Rank) for _ in 1:ClusterCount]
    thread_buffers = [zeros(Float64, length(Stress)) for _ in 1:ThreadCount] # <-- Allocates every call!

    @time StressShear = H2MatVecCalculation(K_d_S, Cluster_W_S, Stress, InAdmissibleBlockSep, ElementRange_SR,  
                                Disp, ThreadCount, thread_buffers,
                                ClusterCount, IsLeaf, IsSingleMom, Cluster_ClusterRange, D_Cluster, 
                                DaughterElem, AdmissibleBlockSep, BlockToClusterSR,  Sig_Cluster, 
                                NonLeafClusterIdx, LeafClusterIdx, Rank, Block_Range_Level)

    Stress = zeros(FaultCount)
    D_Cluster = [zeros(Float64, Rank) for _ in 1:ClusterCount]
    Sig_Cluster = [zeros(Float64, Rank) for _ in 1:ClusterCount]
    thread_buffers = [zeros(Float64, length(Stress)) for _ in 1:ThreadCount] # <-- Allocates every call!
    @time StressNormal = H2MatVecCalculation(K_d_N, Cluster_W_N, Stress, InAdmissibleBlockSep, ElementRange_SR,  
                                Disp, ThreadCount, thread_buffers,
                                ClusterCount, IsLeaf, IsSingleMom, Cluster_ClusterRange, D_Cluster, 
                                DaughterElem, AdmissibleBlockSep, BlockToClusterSR,  Sig_Cluster, 
                                NonLeafClusterIdx, LeafClusterIdx, Rank, Block_Range_Level)














                                
Fig = lines(StressShear)
Fig

##

Fig = lines(StressNormal)
Fig

# GLMakie.activate!()

# # Open window 1
# fig1 = scatter(10.0 .* rand(10, 2))
# display(GLMakie.Screen(), fig1)

# # Open window 2
# fig2 = scatter(10.0 .* rand(10, 2), color=:red)
# display(GLMakie.Screen(), fig2)