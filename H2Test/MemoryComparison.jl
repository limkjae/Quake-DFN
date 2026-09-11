using DelimitedFiles
using Base
using Statistics
using JLD2
using LowRankApprox
using Clustering
using LinearAlgebra
using Base.Threads
using BenchmarkTools

include("../scripts/Functions_H2.jl") 
include("../scripts/Functions_Hmatrix.jl") 

ElementCount = 2537
H2File = "H2Test/Input_Discretized_$(ElementCount)_H2.jld2"
HFile = "H2Test/Input_Discretized_$(ElementCount)_H.jld2"


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

StressShearH2 = zeros(FaultCount)
StressShearF = zeros(FaultCount)
StressShearH = zeros(FaultCount)

StressNormalH2 = zeros(FaultCount)
StressNormalF = zeros(FaultCount)
StressNormalH = zeros(FaultCount)


AdmissibleBlockSep, InAdmissibleBlockSep = 
            Ad_InAd_Blocks_EachThreads(BlockCount, Admissible, ThreadCount, FaultCount, 
                                Block_Range_Level, BlockToClusterSR)


    Stress = zeros(FaultCount)
    D_Cluster = [zeros(Float64, Rank) for _ in 1:ClusterCount]
    Sig_Cluster = [zeros(Float64, Rank) for _ in 1:ClusterCount]
    thread_buffers = [zeros(Float64, length(Stress)) for _ in 1:ThreadCount] # <-- Allocates every call!

    StressShear = H2MatVecCalculation(K_d_S, Cluster_W_S, Stress, InAdmissibleBlockSep, ElementRange_SR,  
                                Disp, ThreadCount, thread_buffers,
                                ClusterCount, IsLeaf, IsSingleMom, Cluster_ClusterRange, D_Cluster, 
                                DaughterElem, AdmissibleBlockSep, BlockToClusterSR,  Sig_Cluster, 
                                NonLeafClusterIdx, LeafClusterIdx, Rank, Block_Range_Level)


                                


    Stress = zeros(FaultCount)
    D_Cluster = [zeros(Float64, Rank) for _ in 1:ClusterCount]
    Sig_Cluster = [zeros(Float64, Rank) for _ in 1:ClusterCount]
    thread_buffers = [zeros(Float64, length(Stress)) for _ in 1:ThreadCount] # <-- Allocates every call!
    StressNormal = H2MatVecCalculation(K_d_N, Cluster_W_N, Stress, InAdmissibleBlockSep, ElementRange_SR,  
                                Disp, ThreadCount, thread_buffers,
                                ClusterCount, IsLeaf, IsSingleMom, Cluster_ClusterRange, D_Cluster, 
                                DaughterElem, AdmissibleBlockSep, BlockToClusterSR,  Sig_Cluster, 
                                NonLeafClusterIdx, LeafClusterIdx, Rank, Block_Range_Level)







################################ H Calculation ############################

FaultCount= load(HFile, "FaultCount")
Ranks_Shear= load(HFile, "Ranks_Shear") # figure(11); plot(Ranks)
Admissible= load(HFile, "Admissible") # figure(11); plot(Ranks)
Ranks_Normal= load(HFile, "Ranks_Normal") # figure(11); plot(Ranks)
ElementRange_SR = load(HFile, "ElementRange_SR")
ShearStiffness_H = load(HFile, "ShearStiffness_H")
NormalStiffness_H = load(HFile, "NormalStiffness_H")
NormalStiffnessZero = load(HFile, "NormalStiffnessZero")



Par_ElementDivision_Shear = ParallelOptimization2(ShearStiffness_H,  
                BlockCount, ThreadCount, Ranks_Shear, Admissible)


StressShearH = HmatSolver_Pararllel(Disp, ShearStiffness_H, ElementRange_SR, 
                                    Par_ElementDivision_Shear, ThreadCount, zeros(FaultCount, ThreadCount))

 StressNormalH = HmatSolver_Pararllel(Disp, NormalStiffness_H, ElementRange_SR, 
                                    Par_ElementDivision_Shear, ThreadCount, zeros(FaultCount, ThreadCount))


println("H-Matrix Done")


println("Shear Stiffness Memory")

Memory_Stiffness_H2 = Base.summarysize(K_d_S) + Base.summarysize(Cluster_W_S)

Memory_Stiffness_H = Base.summarysize(ShearStiffness_H)
# Base.summarysize(K_d_S[55])
# size.(K_d_S)[:]

println(Memory_Stiffness_H2/1e6)
println(Memory_Stiffness_H/1e6)
save("0MemoryH2.jld2", "K_d_S", K_d_S, "Cluster_W_S",Cluster_W_S)
save("0MemoryH.jld2", "ShearStiffness_H",ShearStiffness_H)

println("B Time H2 ")
@btime H2MatVecCalculation(K_d_S, Cluster_W_S, Stress, InAdmissibleBlockSep, ElementRange_SR,  
                            Disp, ThreadCount, thread_buffers,
                            ClusterCount, IsLeaf, IsSingleMom, Cluster_ClusterRange, D_Cluster, 
                            DaughterElem, AdmissibleBlockSep, BlockToClusterSR,  Sig_Cluster, 
                            NonLeafClusterIdx, LeafClusterIdx, Rank, Block_Range_Level)
println("B Time H ")                            
@btime HmatSolver_Pararllel(Disp, ShearStiffness_H, ElementRange_SR, 
                            Par_ElementDivision_Shear, ThreadCount, zeros(FaultCount, ThreadCount))

println("Time Comparison Done")


# Fig = lines(StressShear)
# lines!(vec(StressShearH))
# Fig

##

Fig = lines(StressNormal)
lines!(vec(StressNormalH))
Fig

# GLMakie.activate!()

# # Open window 1
# fig1 = scatter(10.0 .* rand(10, 2))
# display(GLMakie.Screen(), fig1)

# # Open window 2
# fig2 = scatter(10.0 .* rand(10, 2), color=:red)
# display(GLMakie.Screen(), fig2)