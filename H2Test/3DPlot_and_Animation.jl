using DelimitedFiles
using Base
using Statistics
using JLD2
using LowRankApprox
using Clustering
using LinearAlgebra
using Base.Threads
using PyPlot
using PyCall
using JLD2
using LinearAlgebra
using Statistics
pygui(true)

include("../scripts/Functions_Plot.jl")
include("../scripts/Functions_H2.jl") 
include("../scripts/Functions_Hmatrix.jl") 


ElementCount = 2537
H2File = "H2Test/Input_Discretized_$(ElementCount)_H2.jld2"
HFile = "H2Test/Input_Discretized_$(ElementCount)_H.jld2"


FaultCenter, Block_Range_Level, Admissible, ElementRange_SR, Rank, BlockToClusterSR, FaultLengthStrike, FaultLengthDip, FaultStrikeAngle, FaultDipAngle, FaultRakeAngle,
FaultCount, Block_Ctr_Diam, BlockCount, ClusterCount, IsLeaf, IsSingleMom, 
DaughterElem, Cluster_ClusterRange, NonLeafClusterIdx, LeafClusterIdx, SingleMomClusterIdx, RegularClusterIdx, 
Cluster_SkeletonGlobal_S, K_d_S, Cluster_W_S,
Cluster_SkeletonGlobal_N, K_d_N, Cluster_W_N =
load(H2File, "FaultCenter", "Block_Range_Level","Admissible","ElementRange_SR", "Rank","BlockToClusterSR", "FaultLengthStrike", "FaultLengthDip", "FaultStrikeAngle", "FaultDipAngle", "FaultRakeAngle",
    "FaultCount", "Block_Ctr_Diam", "BlockCount", "ClusterCount", "IsLeaf", "IsSingleMom",
    "DaughterElem", "Cluster_ClusterRange", "NonLeafClusterIdx",  "LeafClusterIdx", "SingleMomClusterIdx", "RegularClusterIdx",
    "Cluster_SkeletonGlobal_S","K_d_S",  "Cluster_W_S",  
    "Cluster_SkeletonGlobal_N", "K_d_N", "Cluster_W_N")

println("H2Matrix Loaded")


LoadingFaultCount = 0
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


@time StressShearH = HmatSolver_Pararllel(Disp, ShearStiffness_H, ElementRange_SR, 
                                    Par_ElementDivision_Shear, ThreadCount, zeros(FaultCount, ThreadCount))

@time StressNormalH = HmatSolver_Pararllel(Disp, NormalStiffness_H, ElementRange_SR, 
                                    Par_ElementDivision_Shear, ThreadCount, zeros(FaultCount, ThreadCount))


println("H-Matrix Done")




########################### 3D Plot ############################


ResultName="Result"
FileName="Results/" * ResultName * ".jld2"
FileNameInput="Results/" * ResultName * "_Input.jld2"

#######################################################################################
############################### Figure Configuration ##################################

PlotRotation = [30, -70]
Transparent = 0 # 1 for transparent fault plot. 0 for no-transparency
Edge = 0 # 0 for no element boudary. 1 for plotting element boundary
MinMax_Axis = 0 # 0 for automatically selected axis minimim and maximum 


########################### H2 3D ############################

PlotInput= StressShearH; ColorMinMax = 0
# PlotInput= StressNormal; ColorMinMax = 0

figure(1)
clf()
MaxVaule, MinValue = FaultPlot_3D_Color_General(FaultCenter,FaultLengthStrike, FaultLengthDip,
    FaultStrikeAngle, FaultDipAngle, FaultRakeAngle, PlotInput, 
    PlotRotation, MinMax_Axis, ColorMinMax, Transparent, Edge, LoadingFaultCount)

ax = subplot(projection="3d")
xlabel("x")
ylabel("y")
plotforcbar=  PyPlot.scatter([1,1],[1,1],0.1, [MinValue,MaxVaule], cmap="jet")
colorbar(plotforcbar, pad=0.15)
figure(1).canvas.draw()
    ax.set_aspect("equal")
########################### H2 3D ############################

PlotInput= StressShearH - StressShear; ColorMinMax = 0
# PlotInput= StressNormalH ; ColorMinMax = 0

figure(2)
clf()
MaxVaule, MinValue = FaultPlot_3D_Color_General2(FaultCenter,FaultLengthStrike, FaultLengthDip,
    FaultStrikeAngle, FaultDipAngle, FaultRakeAngle, PlotInput, 
    PlotRotation, MinMax_Axis, ColorMinMax, Transparent, Edge, LoadingFaultCount)

ax = subplot(projection="3d")
xlabel("x")
ylabel("y")
plotforcbar=  PyPlot.scatter([1,1],[1,1],0.1, [MinValue,MaxVaule], cmap="jet")
colorbar(plotforcbar, pad=0.15)
figure(1).canvas.draw()
    ax.set_aspect("equal")

