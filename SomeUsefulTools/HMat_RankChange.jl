
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
using Statistics
@pyimport matplotlib.patches as patches
pygui(true)

############# Only Rank can be changed by adjusting Tolerance #############
##################### Original Matrix Should Exists #######################





######################## Recompression Tolerance ##########################
Tolerance = 1e3 # pascal for 1m slip (More approximaion for higher Tolerance)
PlotHMat = 0 # HMatrix structure plot
###########################################################################

include("../Functions_Solvers.jl")
include("../Functions_RSFDFN3DMain_H.jl")
include("../Results/Functions_Plot.jl")
include("../QuickParameterAdjust.jl")
include("../Functions_Hmatrix.jl")
LoadingInputFileName="Input_Discretized.jld2" 




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
            c = PyObject(patches.Rectangle((ElementRange_SR[i,3]-1, -ElementRange_SR[i,1]+1), ElementRange_SR[i,4] - ElementRange_SR[i,3]+1, 
                        -ElementRange_SR[i,2] + ElementRange_SR[i,1]-1, linewidth=1, edgecolor="k", facecolor=[0.4  0.4  1]))    
        else           
            c = PyObject(patches.Rectangle((ElementRange_SR[i,3]-1, -ElementRange_SR[i,1]+1), ElementRange_SR[i,4] - ElementRange_SR[i,3]+1, 
                        -ElementRange_SR[i,2] + ElementRange_SR[i,1]-1, linewidth=1, edgecolor="k", facecolor=[1 0.4 0.4]))
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
            OrigianlMatrixToApproximate = StiffnessMatrixShear[ElementRange_SR[i,1]:ElementRange_SR[i,2],ElementRange_SR[i,3]:ElementRange_SR[i,4]]
            ApproxMatrixS = pqrfact(OrigianlMatrixToApproximate, atol = Tolerance)
            push!(ShearStiffness_H,ApproxMatrixS)
            Ranks_Shear[BlockIndex] = size(ApproxMatrixS[:Q],2)
            
            OrigianlMatrixToApproximate = StiffnessMatrixNormal[ElementRange_SR[i,1]:ElementRange_SR[i,2],ElementRange_SR[i,3]:ElementRange_SR[i,4]]
            ApproxMatrixN = pqrfact(OrigianlMatrixToApproximate, atol = Tolerance)
            push!(NormalStiffness_H,ApproxMatrixN)
            Ranks_Normal[BlockIndex] = size(ApproxMatrixN[:Q],2)
        else 
            push!(ShearStiffness_H,StiffnessMatrixShear[ElementRange_SR[i,1]:ElementRange_SR[i,2],ElementRange_SR[i,3]:ElementRange_SR[i,4]])
            push!(NormalStiffness_H,StiffnessMatrixNormal[ElementRange_SR[i,1]:ElementRange_SR[i,2],ElementRange_SR[i,3]:ElementRange_SR[i,4]])
        end

    end
    ShearStiffness_H = ShearStiffness_H[2:end]
    NormalStiffness_H = NormalStiffness_H[2:end]

########################## end of HMatrix Approximation ##########################

figure(11)
plot(Ranks_Shear,"r")
figure(12)
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



# # jldsave("HmatSave.jld2"; StiffnessMatrixShearOriginal, StiffnessMatrixNormalOriginal, Ranks, ElementRange_SR)
# Admissible = load(LoadingInputFileName,"Admissible")
# FaultCount= load(LoadingInputFileName, "FaultCount")
# Ranks_Shear= load(LoadingInputFileName, "Ranks_Shear") # figure(11); plot(Ranks)
# Ranks_Normal= load(LoadingInputFileName, "Ranks_Normal") # figure(11); plot(Ranks)
# ElementRange_SR = load(LoadingInputFileName, "ElementRange_SR")
# ShearStiffness_H = load(LoadingInputFileName, "ShearStiffness_H")
# NormalStiffness_H = load(LoadingInputFileName, "NormalStiffness_H")
# NormalStiffnessZero = load(LoadingInputFileName, "NormalStiffnessZero")
# StiffnessMatrixShear= load(LoadingInputFileName, "StiffnessMatrixShear")
# StiffnessMatrixNormal= load(LoadingInputFileName, "StiffnessMatrixNormal")


# ################ Optimizing and Initializing Hmatrix ################
# ThreadCount = Threads.nthreads()
# Epsilon_MaxDiffRatio = 1e-7
# MaxRatioAllowed = 1.5
# MaxIteration = 50

# BlockCount = length(Ranks_Shear)
# println("Shear Stiffness")
# Par_ElementDivision_Shear = ParallelOptimization(ShearStiffness_H, ElementRange_SR, 
#                             FaultCount, BlockCount, ThreadCount, MaxRatioAllowed, MaxIteration)
# println("Normal Stiffness")
# Par_ElementDivision_Normal = ParallelOptimization(NormalStiffness_H, ElementRange_SR, 
#                             FaultCount, BlockCount, ThreadCount, MaxRatioAllowed, MaxIteration)

# #####################################################################

# TestVector = ones(FaultCount)

# for i=1:20
# @time SoluationHMat = HmatSolver_Pararllel(TestVector, ShearStiffness_H, ElementRange_SR, FaultCount,
#                                      Par_ElementDivision_Shear, ThreadCount)
# # println("end")
# end

# for i=1:20
# @time SoluationOrignial = StiffnessMatrixShear * TestVector
# end

# for i=1:5
#     @time SoluationHMat = HmatSolver_ThreadTime(TestVector, ShearStiffness_H, ElementRange_SR, FaultCount,
#                                          Par_ElementDivision_Shear, ThreadCount)
#     println("done")
# end
    
# println("end")
# # plot(SoluationHMat)
# # plot(SoluationOrignial)


# # Ranks= load(LoadingInputFileName, "Ranks")
# # ElementRange_SR = load(LoadingInputFileName, "ElementRange_SR")
# # ShearStiffness_H = load(LoadingInputFileName, "ShearStiffness_H")
# # FaultCount= load(LoadingInputFileName, "FaultCount")
# # BlockCount = length(Ranks)
# # Elastic_Load_DispP = zeros(FaultCount)
# # Elastic_Load_Disp = zeros(FaultCount)
# # NetDisp = ones(FaultCount)
# # # ElasticLoadingShearMatrix = load(LoadingInputFileName, "StiffnessMatrixShear")
# # StiffnessMatrixShear= load(LoadingInputFileName, "StiffnessMatrixShear")


# # BlockNumber = 1
# # AMatrix_Original = StiffnessMatrixShear[ElementRange_SR[BlockNumber,1]:ElementRange_SR[BlockNumber,2],
# #                                         ElementRange_SR[BlockNumber,3]:ElementRange_SR[BlockNumber,4]]

# # AMatrix_QR = ShearStiffness_H[BlockNumber]
# # AMatrix_Sketch = sketchfact(AMatrix_Original, rank=7)
# # AMatrix_Sketch = idfact(AMatrix_Original,  sketch=:srft)
# # test = ones(3958)

# # AMatrix_Original * test
# # AMatrix_QR * test
# # AMatrix_Sketch * test

# # ShearStiffness_H[1][:Q] = ShearStiffness_H[1][:Q][:,1:3]
# # ShearStiffness_H[1][:Q] = 1
# # @time for repeat = 1:50
# # ShearStiffness_H[1][:Q][:,1:3] * ShearStiffness_H[1][:R][1:3,:] * test
# # end
# # @time for repeat = 1:50
# # ShearStiffness_H[1] * test
# # end

# # ArrangePoint = [0,-5000,1000]
# # HowManyDivisionEachLevel = 2
# # TotalHierarchyLevel = 8
# # MinimumElementsToCut = 50
# # ElementPartRoughCount = 1000
# # DistDiamRatioCrit = 1
# # AllowedError = 1e5 # pascal for 1m slip
# # # AllowedError = 1e3 # pascal for 1m slip


# #     ################################## Approximate ######################################
# #     BlockCount = length(ElementRange_SR[:,1])
# #     ShearStiffness_H = Any[0]
# #     BlockIndex = 0
# #     Ranks[Ranks .> 0] .= 1
# #     println("compressing")
# #     for i=1:BlockCount
# #         BlockIndex = BlockIndex + 1
        
# #         if Ranks[BlockIndex] > 0
# #             OrigianlMatrixToApproximate = StiffnessMatrixShear[ElementRange_SR[i,1]:ElementRange_SR[i,2],ElementRange_SR[i,3]:ElementRange_SR[i,4]]
# #             TestVector = ones(size(OrigianlMatrixToApproximate,2))
# #             Error = 1e9
# #             ApproxMatrix = pqrfact(OrigianlMatrixToApproximate, rank=1)
# #             # while Error > 1e2
# #             while Error > AllowedError
# #                 Ranks[BlockIndex] += 1
# #                 ApproxMatrix = pqrfact(OrigianlMatrixToApproximate, rank=Ranks[BlockIndex])
# #                 TestMul_Original = OrigianlMatrixToApproximate * TestVector
# #                 TestMul_Approx = ApproxMatrix * TestVector
# #                 Error = maximum(abs.((TestMul_Original .- TestMul_Approx)))
# #                 # println(Error)
# #             end
# #             push!(ShearStiffness_H,pqrfact(OrigianlMatrixToApproximate, rank=10))
            
# #         else 
# #             push!(ShearStiffness_H,StiffnessMatrixShear[ElementRange_SR[i,1]:ElementRange_SR[i,2],ElementRange_SR[i,3]:ElementRange_SR[i,4]])
# #         end

# #     end
# #     ShearStiffness_H = ShearStiffness_H[2:end]





# # ThreadCount = 24
# # BlockCount = length(Ranks)
# # Par_ElementDivision = ParallelOptimization(ShearStiffness_H, ElementRange_SR, FaultCount, BlockCount, ThreadCount)
# # Elastic_Load_DispP = zeros(FaultCount)
# # ElasticLoadDisp = zeros(FaultCount)

# # @time for i=1:10
# #  Elastic_Load_DispP= 
# #     HmatSolver_Pararllel(NetDisp, ShearStiffness_H, ElementRange_SR, FaultCount, Par_ElementDivision, ThreadCount)
# # end

# # @time for i=1:10
# #  ElasticLoadDisp = StiffnessMatrixShear * NetDisp
# # end
# # maximum(abs.((ElasticLoadDisp - Elastic_Load_DispP)./ElasticLoadDisp))
# # figure(1)
# # plot(ElasticLoadDisp)
# # plot(Elastic_Load_DispP)

# # # @time for i=1:10
# # #     # @time Elastic_Load_Disp  = HmatSolver(NetDisp, ShearStiffness_H, BlockCount, ElementRange_SR, FaultCount)
# # #     @time Elastic_Load_DispP= 
# # #         HmatSolver_Pararllel(NetDisp, ShearStiffness_H, ElementRange_SR, FaultCount, Par_ElementDivision, ThreadCount) ./ Elastic_Load_DispP
# # #     @time ElasticLoadDisp = StiffnessMatrixShear * NetDisp
# # #         # println("One Step Over")
# # # end


# # # Elastic_Load_DispP

# # # println(Par_ElementDivision)
# # # SizeAndRankCumNorm = zeros(size(ShearStiffness_H,1))
# # # SizeAndRank = zeros(size(ShearStiffness_H,1))
# # # for i in eachindex(ShearStiffness_H[:])
# # #     if Ranks[i]>0
# # #         SizeAndRank[i] = size(ShearStiffness_H[i],1) * Ranks[i]
# # #     else 
# # #         SizeAndRank[i] = size(ShearStiffness_H[i],1)^2 
# # #     end

# # #     if i==1; SizeAndRankCumNorm[i] = SizeAndRank[i]
# # #     else
# # #         SizeAndRankCumNorm[i] =  SizeAndRankCumNorm[i-1] + SizeAndRank[i]
# # #     end
# # # end

# # # SizeAndRankCumNorm= SizeAndRankCumNorm/maximum(SizeAndRankCumNorm)


# # # # RecordedTime = ParallelOptimization(ShearStiffness_H, ElementRange_SR, FaultCount, BlockCount)
# # # CumRecordTime = cumsum(RecordedTime)
# # # CumRecordTimeNorm = CumRecordTime / maximum(CumRecordTime)
# # # # plot(CumRecordTimeNorm)

# # # ThreadCount = 24
# # # ElementGap = 1 / ThreadCount
# # # ElementCut = zeros(Int, ThreadCount +1)
# # # for i=1:ThreadCount
# # #     if i==ThreadCount 
# # #         ElementCut[i+1] = BlockCount
# # #     else 
# # #         ElementCut[i+1] = findmin(abs.(CumRecordTimeNorm .- ElementGap * i))[2]
# # #     end
# # # end
# # # Timeone = 0
# # # TimeP = 0

# # # println(ElementCut)
# # # println(Timeone)
# # # println(TimeP)
    
# # # findmin(abs.(SizeAndRankCumNorm .- 0.1))[2]


# # # plot(SizeAndRankCumNorm)
# # # ElementCut = zeros(Int,ThreadCount+1)
# # # for i=1:ThreadCount
# # #     if i==ThreadCount 
# # #         ElementCut[i+1] = BlockCount
# # #     else 
# # #         ElementCut[i+1] = ElementCut[i]  + round(10 * 1.4^i)
# # #     end
# # # end




# # # ElementCut = [0,20,60,90,350,1250,BlockCount]

# # # ElementCut = [0,10,20,35,50,65,85,900,BlockCount]
# # # ElementCut = [0,10,20,35,50,65,85,900,BlockCount]
# # # ElementCut = [0,2,40,60,100,350,550,1550,2800,BlockCount]
# # # ElementCut = [0,20,50,90,150,500,1000,2000,BlockCount]

# # # ThreadCount = length(ElementCut) - 1 

# # # println(maximum(abs.(Elastic_Load_Disp - Elastic_Load_DispP)))

# # # findmin(abs.(ElementCut .- 250))[2]
# # # Elastic_Load_Disp - Elastic_Load_DispP


# # # plot(Elastic_Load_Disp - Elastic_Load_DispP)


# # # Threads.nthreads()
# # # for i=1:10
# # #    println(size(ShearStiffness_H[i],1))
# # # end

# # # SizeDistribution = zeros(size(ShearStiffness_H,1))
# # # for i=1:size(ShearStiffness_H,1)
# # #     if i==1; SizeDistribution[i] = size(ShearStiffness_H[i],1);
# # #     else
# # #         SizeDistribution[i] =  SizeDistribution[i-1] + size(ShearStiffness_H[i],1)
# # #     end
# # # # size(ShearStiffness_H[1],1)

# # # end
# # # plot(SizeDistribution)
# # # Elastic_Load_Disp - Elastic_Load_DispP

# # #  for i=1:2
# # #   
# # # end

# # # @time Elastic_Load_Disp{2, Elastic_Load_Disp2, Elastic_Load_Disp3,Elastic_Load_Disp4} = 
# # #     HmatSolver_Pararllel(NetDisp, ShearStiffness_H, BlockCount, ElementRange_SR, FaultCount)
# # #      Elastic_Load_DispP = Elastic_Load_Disp1 + Elastic_Load_Disp2 + Elastic_Load_Disp3 +Elastic_Load_Disp4

# # # Elastic_Load_DispP - Elastic_Load_Disp
# # # ShearStiffness_H, ElasticLoadingShearMatrix, ReducedStiffnessMatrixShear, ReducedStiffnessMatrixNormal, Ranks, ElementRange_SR = 
# # #     load("HmatSave.jld2", "ShearStiffness_H", "ElasticLoadingShearMatrix", "ReducedStiffnessMatrixShear", "ReducedStiffnessMatrixNormal",
# # #          "Ranks", "ElementRange_SR")
# # # # FaultCount = length(StiffnessMatrixShear[:,1])

# # # NetDisp = ones(FaultCount)
# # # Elastic_Load_Disp = zeros(FaultCount)
# # # for Blockidx in 1:BlockCount             
# # #     Elastic_Load_Disp[ElementRange_SR[Blockidx,1]:ElementRange_SR[Blockidx,2]] = 
# # #         Elastic_Load_Disp[ElementRange_SR[Blockidx,1]:ElementRange_SR[Blockidx,2]] + 
# # #         ShearStiffness_H[Blockidx] * NetDisp[ElementRange_SR[Blockidx,3]:ElementRange_SR[Blockidx,4]]
# # # end

# # # @time Elastic_Load_DispOriginal=ElasticLoadingShearMatrix * (NetDisp) 

# # # Elastic_Load_DispOriginal-Elastic_Load_Disp

# # # # Ranks = RanksOriginal

