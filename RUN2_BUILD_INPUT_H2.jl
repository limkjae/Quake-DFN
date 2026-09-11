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
# @pyimport matplotlib.patches as patches
pygui(true)

# include("FindAwayElements.jl")
# include("Functions_Plot.jl")


include("scripts/Functions_H2.jl") 
include("scripts/Functions_BuildInputFile.jl")
include("scripts/Functions_OKADA3D.jl")
# include("scripts/Functions_OKADA3D_Vec.jl")
include("scripts/Functions_Plot.jl")
include("scripts/Functions_Hmatrix.jl")


########################## input ############################
InputBulkFileName="Input_BulkFaultGeometry.txt"
OutputFileName="Input_Discretized.jld2"
HMatrixStructureFile = "Input_HmatrixStructure.jld2"

Rank = 20
RadMultipleForElemSelect = 4.0
MinimumProxyElments = 200
if MinimumProxyElments < Rank; MinimumProxyElments = Rank; end
#############################################################


########################### Load the Input Data #############################

Block_Ctr_Diam, Block_Range_Level, Input_Segment, LoadingFaultExist, LoadingFaultCount,
MinimumElementsToCut, ArrangePoint, Admissible, ElementRange_SR, Switch_StrikeSlip_or_ReverseNormal, 
ShearModulus, PoissonRatio, RockDensity, DropCrit, DropCritNormalStressMultiplier, MinimumNS, RorT, DistDiamRatioCrit =
load(HMatrixStructureFile, "Block_Ctr_Diam", "Block_Range_Level", "Input_Segment", "LoadingFaultExist", "LoadingFaultCount",
    "MinimumElementsToCut", "ArrangePoint", "Admissible", "ElementRange_SR", "Switch_StrikeSlip_or_ReverseNormal", 
    "ShearModulus", "PoissonRatio", "RockDensity", "DropCrit", "DropCritNormalStressMultiplier", "MinimumNS", "RorT","DistDiamRatioCrit")


FaultCount = length(Input_Segment[:,1])

######################## Read Segmented Input File #######################
FaultCenter, Fault_a, Fault_b, Fault_Dc, Fault_Theta_i, Fault_V_i, Fault_Friction_i, Fault_NormalStress, 
Fault_V_Const, Fault_BulkIndex, FaultLengthStrike, FaultLengthDip, FaultStrikeAngle, 
FaultDipAngle, FaultRakeAngle, FaultLengthStrike_Bulk, FaultLengthDip_Bulk, NormalStiffnessZero = 
    ReadSegmentInput(Input_Segment, FaultCount, RorT) 



############ Find the source and receiver clusters for each block ###########

Block_Range_Level = Block_Range_Level[2:end,:] # Remove the loading elements
Block_Ctr_Diam = Block_Ctr_Diam[2:end,:] # Remove the loading elements
Reorder = sortperm(Block_Range_Level[:,3], rev=true) # Move the low-level (large) block toward top
Block_Range_Level = Block_Range_Level[Reorder,:]
Block_Ctr_Diam = Block_Ctr_Diam[Reorder,:]
ClusterCount = length(Block_Range_Level[:,1])
BlockCount= length(ElementRange_SR[:,1])
FinalLevel = maximum(Block_Range_Level[:,3])

BlockToClusterSR = zeros(Int,BlockCount,2)
ClusterCount = length(Block_Range_Level[:,1])
for BlockIdx = 1:BlockCount
    for ClusterIdx = 1:ClusterCount
        if ElementRange_SR[BlockIdx,1:2] == Block_Range_Level[ClusterIdx,1:2]
            BlockToClusterSR[BlockIdx,1] = ClusterIdx
        end
        if ElementRange_SR[BlockIdx,3:4] == Block_Range_Level[ClusterIdx,1:2]
            BlockToClusterSR[BlockIdx,2] = ClusterIdx
        end
    end
end

######################## Define Mom-Daughter Relationship ##################

 LeafClusterIdx, SingleMomClusterIdx, RegularClusterIdx, NonLeafClusterIdx, IsSingleMom, IsLeaf, DaughterElemVec, DaughterElem = 
    FindMomDaughter(Block_Range_Level, ClusterCount)


####################### Calculate the Skeletons and W  for Leaf #######################


Cluster_SkeletonLocal_S, Cluster_SkeletonGlobal_S, Cluster_W_S, Cluster_ClusterRange,
    Cluster_SkeletonLocal_N, Cluster_SkeletonGlobal_N, Cluster_W_N = 
        FindW_SN(ClusterCount, Input_Segment, Block_Range_Level, DaughterElemVec, IsLeaf, IsSingleMom,
                 RadMultipleForElemSelect, MinimumProxyElments, DistDiamRatioCrit, Rank, ShearModulus, PoissonRatio) 

println("SN Done")



K_d_S, K_d_N = FindKd(BlockCount, ElementRange_SR, Admissible, BlockToClusterSR,
                Cluster_SkeletonLocal_S, Cluster_W_S, Cluster_ClusterRange,
                Cluster_SkeletonLocal_N, Cluster_W_N )




save(OutputFileName, 
"FaultCenter", FaultCenter,
"ShearModulus", ShearModulus, "RockDensity", RockDensity, "PoissonRatio", PoissonRatio,
"FaultLengthStrike", FaultLengthStrike, "FaultLengthDip", FaultLengthDip, "FaultStrikeAngle", FaultStrikeAngle, 
"FaultDipAngle", FaultDipAngle, "FaultRakeAngle", FaultRakeAngle, "Fault_a", Fault_a, "Fault_b", Fault_b, "Fault_Dc", Fault_Dc, 
"Fault_Theta_i", Fault_Theta_i, "Fault_V_i", Fault_V_i, "Fault_Friction_i", Fault_Friction_i, "Fault_NormalStress", Fault_NormalStress, 
"Fault_V_Const", Fault_V_Const, "Fault_BulkIndex", Fault_BulkIndex, "FaultLengthStrike_Bulk", FaultLengthStrike_Bulk, 
"FaultLengthDip_Bulk", FaultLengthDip_Bulk, "FaultCount", FaultCount, "LoadingFaultCount", LoadingFaultCount, 
"Switch_StrikeSlip_or_ReverseNormal", Switch_StrikeSlip_or_ReverseNormal, "MinimumNormalStress", MinimumNS,
"NormalStiffnessZero", NormalStiffnessZero, "RorT", RorT,   
"Block_Ctr_Diam", Block_Ctr_Diam, "Block_Range_Level", Block_Range_Level, "ArrangePoint", ArrangePoint, 
"Admissible", Admissible, "ElementRange_SR", ElementRange_SR, 
  "Rank",Rank, "BlockToClusterSR",BlockToClusterSR,
"BlockCount", BlockCount, "IsLeaf", IsLeaf, "IsSingleMom", IsSingleMom, "DaughterElem", DaughterElem, "ClusterCount",ClusterCount,
"NonLeafClusterIdx", NonLeafClusterIdx, "LeafClusterIdx", LeafClusterIdx, "SingleMomClusterIdx", SingleMomClusterIdx, "RegularClusterIdx", RegularClusterIdx,
"Cluster_ClusterRange", Cluster_ClusterRange, 
"Cluster_SkeletonGlobal_S",Cluster_SkeletonGlobal_S, "K_d_S", K_d_S, "Cluster_W_S", Cluster_W_S, 
"Cluster_SkeletonGlobal_N",Cluster_SkeletonGlobal_N, "K_d_N", K_d_N, "Cluster_W_N", Cluster_W_N)
println("Saved File Name: ",OutputFileName)

