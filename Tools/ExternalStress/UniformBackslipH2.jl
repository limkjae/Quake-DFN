
using DelimitedFiles
using JLD2
using LinearAlgebra
using Printf
using SpecialFunctions
using PyPlot
using PyCall
using LowRankApprox
using Statistics
pygui(true)


LoadingRate=1e-9

OutputFile="Input_ExternalStressChange.jld2"
LoadingInputFileName="Input_Discretized.jld2" 

include("../../scripts/Functions_H2.jl")



FaultCenter, Block_Range_Level, Admissible, ElementRange_SR, Rank, BlockToClusterSR,
FaultCount, Block_Ctr_Diam, BlockCount, ClusterCount, IsLeaf, IsSingleMom, 
DaughterElem, Cluster_ClusterRange, NonLeafClusterIdx, LeafClusterIdx, SingleMomClusterIdx, RegularClusterIdx, 
Cluster_SkeletonGlobal_S, K_d_S, Cluster_W_S,
Cluster_SkeletonGlobal_N, K_d_N, Cluster_W_N =
load(LoadingInputFileName, "FaultCenter", "Block_Range_Level","Admissible","ElementRange_SR", "Rank","BlockToClusterSR",
    "FaultCount", "Block_Ctr_Diam", "BlockCount", "ClusterCount", "IsLeaf", "IsSingleMom",
    "DaughterElem", "Cluster_ClusterRange", "NonLeafClusterIdx",  "LeafClusterIdx", "SingleMomClusterIdx", "RegularClusterIdx",
    "Cluster_SkeletonGlobal_S","K_d_S",  "Cluster_W_S",  
    "Cluster_SkeletonGlobal_N", "K_d_N", "Cluster_W_N")


ThreadCount = 5

Disp = ones(FaultCount) 

StreeVectorShear = zeros(FaultCount)
StreeVectorNormal = zeros(FaultCount)

Stress = zeros(FaultCount)
D_Cluster = [zeros(Float64, Rank) for _ in 1:ClusterCount]
Sig_Cluster = [zeros(Float64, Rank) for _ in 1:ClusterCount]
thread_buffers = [zeros(Float64, length(Stress)) for _ in 1:ThreadCount] # <-- Allocates every call!

StreeVectorShear = H2MatVecCalculation(K_d_S, Cluster_W_S, Stress, InAdmissibleBlockSep, ElementRange_SR,  
                            Disp, ThreadCount, thread_buffers,
                            ClusterCount, IsLeaf, IsSingleMom, Cluster_ClusterRange, D_Cluster, 
                            DaughterElem, AdmissibleBlockSep, BlockToClusterSR,  Sig_Cluster, 
                            NonLeafClusterIdx, LeafClusterIdx, Rank, Block_Range_Level)

Stress = zeros(FaultCount)
D_Cluster = [zeros(Float64, Rank) for _ in 1:ClusterCount]
Sig_Cluster = [zeros(Float64, Rank) for _ in 1:ClusterCount]
thread_buffers = [zeros(Float64, length(Stress)) for _ in 1:ThreadCount] # <-- Allocates every call!
StreeVectorNormal = H2MatVecCalculation(K_d_N, Cluster_W_N, Stress, InAdmissibleBlockSep, ElementRange_SR,  
                            Disp, ThreadCount, thread_buffers,
                            ClusterCount, IsLeaf, IsSingleMom, Cluster_ClusterRange, D_Cluster, 
                            DaughterElem, AdmissibleBlockSep, BlockToClusterSR,  Sig_Cluster, 
                            NonLeafClusterIdx, LeafClusterIdx, Rank, Block_Range_Level)

ExternalStress_Normal = zeros(2,FaultCount)
ExternalStress_Shear =  zeros(2,FaultCount)
ExternalStress_Shear[2,:] = - StreeVectorShear* 1e11 * 1e-9
ExternalStress_Normal[2,:] =  - StreeVectorNormal* 1e11 * 1e-9
ExternalStress_TimeArray = zeros(2)
ExternalStress_TimeArray[2] = 1e11
Pressure =  zeros(2,FaultCount)
FlowRate = 0
PressureOrigin = [0,0,0]


save(OutputFile, 
"ExternalStress_Normal", ExternalStress_Normal, "ExternalStress_Shear", ExternalStress_Shear, "Pressure", Pressure,
"ExternalStress_TimeArray", ExternalStress_TimeArray, "FlowRate", FlowRate, "PressureOrigin", PressureOrigin)
    
