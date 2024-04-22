
using PyPlot
using PyCall
using DelimitedFiles
using JLD2
using LinearAlgebra
using Printf
using SpecialFunctions
using HMatrices
using StaticArrays
using LowRankApprox
using Distributed
using LoopVectorization
pygui(true)


include("Functions_Solvers.jl")
include("Functions_RSFDFN3DMain_H.jl")
include("Results/Functions_Plot.jl")
include("QuickParameterAdjust.jl")
include("Functions_Hmatrix.jl")
LoadingInputFileName="Input_Discretized_H.jld2" 

# jldsave("HmatSave.jld2"; StiffnessMatrixShearOriginal, StiffnessMatrixNormalOriginal, Ranks, ElementRange_SR)

StiffnessMatrixShear= load(LoadingInputFileName, "StiffnessMatrixShear")
StiffnessMatrixNormal= load(LoadingInputFileName, "StiffnessMatrixNormal")
Ranks= load(LoadingInputFileName, "Ranks")
ElementRange_SR = load(LoadingInputFileName, "ElementRange_SR")
ShearStiffness_H = load(LoadingInputFileName, "ShearStiffness_H")
FaultCount = size(StiffnessMatrixShear,1)

InitialNormalStress= load(LoadingInputFileName, "Fault_NormalStress")
FrictionI= load(LoadingInputFileName, "Fault_Friction_i")
InitialShearStress = InitialNormalStress .* FrictionI

BlockCount = length(Ranks)
Elastic_Load_DispP = zeros(FaultCount)
Elastic_Load_Disp = zeros(FaultCount)
NetDisp = ones(FaultCount)
FaultCount= size(StiffnessMatrixShear,1)

# K_Self=zeros(FaultCount)
# for i=1:FaultCount
#     K_Self[i]=abs(StiffnessMatrixShear[i,i])
# end
# Far_Load_Disp_Initial_CorrectSolution = -(StiffnessMatrixShear\InitialShearStress)
# StiffnessMatrixShear * Far_Load_Disp_Initial_CorrectSolution

Epsilon_MaxDiffRatio = 1e-5

ThreadCount = 20
Par_ElementDivision = ParallelOptimization(ShearStiffness_H, ElementRange_SR, FaultCount, BlockCount, ThreadCount)


LoadingStiffnessH, K_Self = StiffnessTransitionToLoading(ShearStiffness_H, ElementRange_SR, FaultCount)


InitialDisplacemenet = SolveAx_b(LoadingStiffnessH, K_Self,  InitialShearStress, ElementRange_SR, FaultCount, 
            Par_ElementDivision, ThreadCount, Epsilon_MaxDiffRatio)

InitialStress = -HmatSolver_Pararllel(InitialDisplacemenet, ShearStiffness_H, ElementRange_SR, FaultCount, Par_ElementDivision, ThreadCount)

InitialShearStress

plot(Far_Load_Disp_Initial_CorrectSolution)
plot(InitialDisplacemenet)



# DispI_k = ones(FaultCount)
# MaxDiff = 1
# iteration=0
# while MaxDiff > 1e-6
#     iteration = iteration+1
#     println(MaxDiff)
# EFTerm = HmatSolver_Pararllel(DispI_k, LoadingStiffnessH, ElementRange_SR, FaultCount, Par_ElementDivision, ThreadCount)

# DispI_kp1 = (EFTerm + InitialShearStress) ./ K_Self
# MaxDiff = maximum(abs.((DispI_k - DispI_kp1) ./ DispI_k))
# DispI_k=copy(DispI_kp1)
# end

# plot(Far_Load_Disp_Initial_CorrectSolution)
# plot(DispI_k)


# # plot(Far_Load_Disp_Initial_CorrectSolution)


# # TestV = ones(FaultCount)


# # Elastic_Load_DispP= 
# #     HmatSolver_Pararllel(-DispI_k, ShearStiffness_H, ElementRange_SR, FaultCount, Par_ElementDivision, ThreadCount)
# #     InitialShearStress
# #     Elastic_Load_DispP-InitialShearStress
# # Elastic_Load_Disp = StiffnessMatrixShear * DispI_k



# #=
 

# UpdagedGuess = zeros(FaultCount)

# Ds= 1e-10
# InitialGuess = Far_Load_Disp_Initial_CorrectSolution .-0.1


# for iteration =1:10
# println(iteration)
# Elastic_Load_DispP_Ini = 
#         HmatSolver_Pararllel(InitialGuess, ShearStiffness_H, ElementRange_SR, FaultCount, Par_ElementDivision, ThreadCount)
# StressDiffIni = (Elastic_Load_DispP_Ini - InitialShearStress)

#     for ChangeElemIndex = 1:FaultCount

#     TestGuess = copy(InitialGuess)
#     TestGuess[ChangeElemIndex] = InitialGuess[ChangeElemIndex] +  Ds

#     Elastic_Load_DispP_After = 
#             HmatSolver_Pararllel(TestGuess, ShearStiffness_H, ElementRange_SR, FaultCount, Par_ElementDivision, ThreadCount)
#     StressDiffAfter = (Elastic_Load_DispP_After - InitialShearStress)

#     Slope = norm(StressDiffIni - StressDiffAfter) / Ds

#     UpdagedGuess[ChangeElemIndex] = InitialGuess[ChangeElemIndex] + StressDiffIni[ChangeElemIndex] / Slope
#     end
# # Far_Load_Disp_Initial_CorrectSolution - UpdagedGuess
# InitialGuess = copy(UpdagedGuess)
# println(maximum(abs.(Far_Load_Disp_Initial_CorrectSolution - UpdagedGuess)))
# end
# #=


# Elastic_Load_DispP_Ini= 
#         HmatSolver_Pararllel(InitialGuess, ShearStiffness_H, ElementRange_SR, FaultCount, Par_ElementDivision, ThreadCount)
# StressDiffInit = Elastic_Load_DispP_Ini - InitialShearStress
# InitialGuessD = InitialGuess .+ Ds
# Elastic_Load_DispP_After= 
#         HmatSolver_Pararllel(InitialGuessD, ShearStiffness_H, ElementRange_SR, FaultCount, Par_ElementDivision, ThreadCount)
# StressDiffAfter = Elastic_Load_DispP_After - InitialShearStress
# Diff = StressDiffInit - StressDiffAfter
# Slope = Diff./Ds
# InitialGuess=InitialGuess .+   StressDiffInit ./Slope

# TestVector = ones(FaultCount)

# Inversed = StiffnessMatrixShear \ TestVector
# StiffnessMatrixShear * Inversed
# ShearStiffness_H[1]
# ################## Numerical Inversion ##################

# InitialGuess = 

# ThreadCount = 20
# Par_ElementDivision = ParallelOptimization(ShearStiffness_H, ElementRange_SR, FaultCount, BlockCount, ThreadCount)
# Elastic_Load_DispP= 
#         HmatSolver_Pararllel(Inversed, ShearStiffness_H, ElementRange_SR, FaultCount, Par_ElementDivision, ThreadCount)




# # TestInv = ones(size(ShearStiffness_H[1],1))
# # TestInv \ ShearStiffness_H[1]
# # size(ShearStiffness_H[1])
# # inv(ShearStiffness_H[1])

# LoadingStiffnessH = StiffnessTransitionToLoading(ShearStiffness_H, ElementRange_SR)

# K_Self=zeros(FaultCount)
# ElasticLoadingShearMatrix=zeros(FaultCount,FaultCount)
# for i=1:FaultCount
#     K_Self[i]=abs(StiffnessMatrixShear[i,i]);
#     for j=1:FaultCount  
#         if i==j
#             ElasticLoadingShearMatrix[i,j]=0.0;
#         else
#             # ElasticLoadingShearMatrix[i,j]=StiffnessMatrixShear[i,j]
#             ElasticLoadingShearMatrix[i,j]=StiffnessMatrixShear[i,j]/StiffnessMatrixShear[i,i];
            
#         end
#     end
# end



# ThreadCount = 20
# Par_ElementDivision = ParallelOptimization(ShearStiffness_H, ElementRange_SR, FaultCount, BlockCount, ThreadCount)
# Elastic_Load_DispP= 
#         HmatSolver_Pararllel(TestVector, ShearStiffness_H, ElementRange_SR, FaultCount, Par_ElementDivision, ThreadCount)

# Elastic_Load_DispP = Elastic_Load_DispP ./ -K_Self

# OriginalSolution = ElasticLoadingShearMatrix * TestVector
        
# Elastic_Load_DispP - OriginalSolution

# #=

# ThreadCount = 20

# BlockCount = length(Ranks)
# Par_ElementDivision = ParallelOptimization(ShearStiffness_H, ElementRange_SR, FaultCount, BlockCount, ThreadCount)


# @time for i=1:30
#     # @time Elastic_Load_Disp  = HmatSolver(NetDisp, ShearStiffness_H, BlockCount, ElementRange_SR, FaultCount)
#     Elastic_Load_DispP= 
#         HmatSolver_Pararllel(NetDisp, ShearStiffness_H, ElementRange_SR, FaultCount, Par_ElementDivision, ThreadCount)
#         #  Elastic_Load_DispP = Elastic_Load_Disp1 + Elastic_Load_Disp2 + Elastic_Load_Disp3 +Elastic_Load_Disp4
#         # println("One Step Over")
# end


# # Elastic_Load_DispP

# # println(Par_ElementDivision)
# # SizeAndRankCumNorm = zeros(size(ShearStiffness_H,1))
# # SizeAndRank = zeros(size(ShearStiffness_H,1))
# # for i in eachindex(ShearStiffness_H[:])
# #     if Ranks[i]>0
# #         SizeAndRank[i] = size(ShearStiffness_H[i],1) * Ranks[i]
# #     else 
# #         SizeAndRank[i] = size(ShearStiffness_H[i],1)^2 
# #     end

# #     if i==1; SizeAndRankCumNorm[i] = SizeAndRank[i]
# #     else
# #         SizeAndRankCumNorm[i] =  SizeAndRankCumNorm[i-1] + SizeAndRank[i]
# #     end
# # end

# # SizeAndRankCumNorm= SizeAndRankCumNorm/maximum(SizeAndRankCumNorm)


# # # RecordedTime = ParallelOptimization(ShearStiffness_H, ElementRange_SR, FaultCount, BlockCount)
# # CumRecordTime = cumsum(RecordedTime)
# # CumRecordTimeNorm = CumRecordTime / maximum(CumRecordTime)
# # # plot(CumRecordTimeNorm)

# # ThreadCount = 24
# # ElementGap = 1 / ThreadCount
# # ElementCut = zeros(Int, ThreadCount +1)
# # for i=1:ThreadCount
# #     if i==ThreadCount 
# #         ElementCut[i+1] = BlockCount
# #     else 
# #         ElementCut[i+1] = findmin(abs.(CumRecordTimeNorm .- ElementGap * i))[2]
# #     end
# # end
# # Timeone = 0
# # TimeP = 0

# # println(ElementCut)
# # println(Timeone)
# # println(TimeP)
    
# # findmin(abs.(SizeAndRankCumNorm .- 0.1))[2]


# # plot(SizeAndRankCumNorm)
# # ElementCut = zeros(Int,ThreadCount+1)
# # for i=1:ThreadCount
# #     if i==ThreadCount 
# #         ElementCut[i+1] = BlockCount
# #     else 
# #         ElementCut[i+1] = ElementCut[i]  + round(10 * 1.4^i)
# #     end
# # end




# # ElementCut = [0,20,60,90,350,1250,BlockCount]

# # ElementCut = [0,10,20,35,50,65,85,900,BlockCount]
# # ElementCut = [0,10,20,35,50,65,85,900,BlockCount]
# # ElementCut = [0,2,40,60,100,350,550,1550,2800,BlockCount]
# # ElementCut = [0,20,50,90,150,500,1000,2000,BlockCount]

# # ThreadCount = length(ElementCut) - 1 

# # println(maximum(abs.(Elastic_Load_Disp - Elastic_Load_DispP)))

# # findmin(abs.(ElementCut .- 250))[2]
# # Elastic_Load_Disp - Elastic_Load_DispP


# # plot(Elastic_Load_Disp - Elastic_Load_DispP)


# # Threads.nthreads()
# # for i=1:10
# #    println(size(ShearStiffness_H[i],1))
# # end

# # SizeDistribution = zeros(size(ShearStiffness_H,1))
# # for i=1:size(ShearStiffness_H,1)
# #     if i==1; SizeDistribution[i] = size(ShearStiffness_H[i],1);
# #     else
# #         SizeDistribution[i] =  SizeDistribution[i-1] + size(ShearStiffness_H[i],1)
# #     end
# # # size(ShearStiffness_H[1],1)

# # end
# # plot(SizeDistribution)
# # Elastic_Load_Disp - Elastic_Load_DispP

# #  for i=1:2
# #   
# # end

# # @time Elastic_Load_Disp{2, Elastic_Load_Disp2, Elastic_Load_Disp3,Elastic_Load_Disp4} = 
# #     HmatSolver_Pararllel(NetDisp, ShearStiffness_H, BlockCount, ElementRange_SR, FaultCount)
# #      Elastic_Load_DispP = Elastic_Load_Disp1 + Elastic_Load_Disp2 + Elastic_Load_Disp3 +Elastic_Load_Disp4

# # Elastic_Load_DispP - Elastic_Load_Disp
# # ShearStiffness_H, ElasticLoadingShearMatrix, ReducedStiffnessMatrixShear, ReducedStiffnessMatrixNormal, Ranks, ElementRange_SR = 
# #     load("HmatSave.jld2", "ShearStiffness_H", "ElasticLoadingShearMatrix", "ReducedStiffnessMatrixShear", "ReducedStiffnessMatrixNormal",
# #          "Ranks", "ElementRange_SR")
# # # FaultCount = length(StiffnessMatrixShear[:,1])

# # NetDisp = ones(FaultCount)
# # Elastic_Load_Disp = zeros(FaultCount)
# # for Blockidx in 1:BlockCount             
# #     Elastic_Load_Disp[ElementRange_SR[Blockidx,1]:ElementRange_SR[Blockidx,2]] = 
# #         Elastic_Load_Disp[ElementRange_SR[Blockidx,1]:ElementRange_SR[Blockidx,2]] + 
# #         ShearStiffness_H[Blockidx] * NetDisp[ElementRange_SR[Blockidx,3]:ElementRange_SR[Blockidx,4]]
# # end

# # @time Elastic_Load_DispOriginal=ElasticLoadingShearMatrix * (NetDisp) 

# # Elastic_Load_DispOriginal-Elastic_Load_Disp

# # # Ranks = RanksOriginal

# =#
# =#
# =#