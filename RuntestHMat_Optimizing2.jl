
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
pygui(true)


include("Functions_Solvers.jl")
include("Functions_RSFDFN3DMain_H.jl")
include("Results/Functions_Plot.jl")
include("QuickParameterAdjust.jl")
include("Functions_Hmatrix.jl")
LoadingInputFileName="Input_Discretized_H_16k.jld2" 

# jldsave("HmatSave.jld2"; StiffnessMatrixShearOriginal, StiffnessMatrixNormalOriginal, Ranks, ElementRange_SR)

Ranks= load(LoadingInputFileName, "Ranks")
ElementRange_SR = load(LoadingInputFileName, "ElementRange_SR")
ShearStiffness_H = load(LoadingInputFileName, "ShearStiffness_H")
FaultCount= load(LoadingInputFileName, "FaultCount")
BlockCount = length(Ranks)
Elastic_Load_DispP = zeros(FaultCount)
Elastic_Load_Disp = zeros(FaultCount)
NetDisp = ones(FaultCount)
# ElasticLoadingShearMatrix = load(LoadingInputFileName, "StiffnessMatrixShear")
StiffnessMatrixShear= load(LoadingInputFileName, "StiffnessMatrixShear")

ThreadCount = 24
BlockCount = length(Ranks)
MaxRatioAllowed = 1.5
MaxIteration = 50
Par_ElementDivision = ParallelOptimization(ShearStiffness_H, ElementRange_SR, 
                FaultCount, BlockCount, ThreadCount, MaxRatioAllowed, MaxIteration)

ElementCountPerDivision = zeros(Int,ThreadCount)

for i=1:ThreadCount
ElementCountPerDivision[i] = Par_ElementDivision[i+1] - Par_ElementDivision[i] 
end
# sum(ElementCountPerDivision)

figure(1)
clf()
figure(2)
clf()

@time ElapsedTime = HmatSolver_SpeedTest(NetDisp, ShearStiffness_H, ElementRange_SR, FaultCount, Par_ElementDivision, ThreadCount)

Maxindex = argmax(ElapsedTime, dims=1)[1]
AverageTime = mean(ElapsedTime)
figure(1)
plot(ElapsedTime)
plot([0,ThreadCount],[AverageTime,AverageTime])
figure(2)
plot(log10.(ElementCountPerDivision))


Minimum_Par_ElementDivision = zeros(length(Par_ElementDivision))
Minimum_ElapsedTimeRatio = 1e5
iteration=0
Termination = 1
MaxIteration = 50
while Termination ==1
# for iteration =1:20
    iteration=iteration+1
    # println(iteration)
    for i=1:ThreadCount
        # if ElapsedTime[i] > AverageTime * 1.2 || ElapsedTime[i] < AverageTime * 0.8
        if ElapsedTime[i] > AverageTime
            ElapsedOverAverage = ElapsedTime[i] / AverageTime
            ElementCountPerDivision[i] = round(Int,ElementCountPerDivision[i] * (1 - (1- 1/ElapsedOverAverage)*0.5))
        else 
            ElapsedOverAverage = ElapsedTime[i] / AverageTime
            if ElementCountPerDivision[i] ==1 
                ElementCountPerDivision[i] = round(Int,ElementCountPerDivision[i] * (1 - (1- 1/ElapsedOverAverage)))
            else
            ElementCountPerDivision[i] = round(Int,ElementCountPerDivision[i] * (1 - (1- 1/ElapsedOverAverage)*0.1))
            end
        end

        # end
        # if ElapsedTime[i] > AverageTime * 1.3
        #     ElementCountPerDivision[i] = round(Int,ElementCountPerDivision[i] * 0.8)
        # end
        
        # if ElapsedTime[i] < AverageTime * 0.5
        #     ElementCountPerDivision[i] = round(Int,ElementCountPerDivision[i] * 1.5)
        # end
        if ElementCountPerDivision[i] < 1
            ElementCountPerDivision[i] = 1
        end
    end

    # ElementCountPerDivision[Maxindex] = round(Int,ElementCountPerDivision[Maxindex] * 0.8)
    ElementCountPerDivision = round.(Int, ElementCountPerDivision * BlockCount / sum(ElementCountPerDivision))

    for i=1:ThreadCount -1 
        Par_ElementDivision[i+1] = ElementCountPerDivision[i] + Par_ElementDivision[i]
    end
    Par_ElementDivision[end]= BlockCount

    for i=1:ThreadCount
        ElementCountPerDivision[i] = Par_ElementDivision[i+1] - Par_ElementDivision[i] 
    end


    ElapsedTime = HmatSolver_SpeedTest(NetDisp, ShearStiffness_H, ElementRange_SR, FaultCount, Par_ElementDivision, ThreadCount)

    Maxindex = argmax(ElapsedTime, dims=1)[1]
    AverageTime = mean(ElapsedTime)
    maxAverageTimeRatio = maximum(ElapsedTime/AverageTime)
    println(maxAverageTimeRatio)

    if maxAverageTimeRatio < Minimum_ElapsedTimeRatio
        Minimum_ElapsedTimeRatio = copy(maxAverageTimeRatio)
        Minimum_Par_ElementDivision = copy(Par_ElementDivision)
        # println(Minimum_Par_ElementDivision)
    end


    if maxAverageTimeRatio < 1.5
        Termination = 0
        println(Par_ElementDivision)
    elseif iteration == MaxIteration
        Termination = 0
        Par_ElementDivision = copy(Minimum_Par_ElementDivision)
        println(Par_ElementDivision)
    end

    figure(1)
    plot(ElapsedTime)
    plot([0,ThreadCount],[AverageTime,AverageTime])
    figure(2)
    plot(log10.(ElementCountPerDivision))
    
end


@time ElapsedTime = HmatSolver_SpeedTest(NetDisp, ShearStiffness_H, ElementRange_SR, FaultCount, Par_ElementDivision, ThreadCount)






# Par_ElementDivision[7] = 50
# Par_ElementDivision[8] = 100
# Par_ElementDivision[9] = 400
# @time for i=1:100
 @time Elastic_Load_DispP= 
    HmatSolver_Pararllel(NetDisp, ShearStiffness_H, ElementRange_SR, FaultCount, Par_ElementDivision, ThreadCount)
# end
println("done")


# @time for i=1:100
 @time ElasticLoadDisp = StiffnessMatrixShear * NetDisp
# end
maximum(abs.((ElasticLoadDisp - Elastic_Load_DispP)))

#=

####################### Recompress ##########################
AllowedError = 1e3

    ################################## Approximate ######################################
    BlockCount = length(ElementRange_SR[:,1])
    ShearStiffness_H2 = Any[0]
    BlockIndex = 0
    # Ranks = Ranks .* 2
    # Ranks= Ranks
    println("compressing")
    for i=1:BlockCount
        # BlockIndex = BlockIndex + 1
        BlockIndex = i
        println(BlockIndex)
        if Ranks[BlockIndex] > 0
            OrigianlMatrixToApproximate = StiffnessMatrixShear[ElementRange_SR[i,1]:ElementRange_SR[i,2],ElementRange_SR[i,3]:ElementRange_SR[i,4]]
            # TestVector = ones(size(OrigianlMatrixToApproximate,2))
            # Error = 1e9
            # Ranks[BlockIndex] = 1
            # # while Error > 1e2
            # while Error > AllowedError
            #     Ranks[BlockIndex] += 1

            #     # ApproxMatrix = pqrfact(OrigianlMatrixToApproximate, rank=Ranks[BlockIndex])
            #     ApproxMatrix = pqrfact(OrigianlMatrixToApproximate, atol = 1e3)
            #     TestMul_Original = OrigianlMatrixToApproximate * TestVector
            #     TestMul_Approx = ApproxMatrix * TestVector
            #     Error = maximum(abs.(TestMul_Original-TestMul_Approx))
            #     # Error2 = maximum(abs.((TestMul_Original .- TestMul_Approx)))
            #     # Error = maximum(abs.(OrigianlMatrixToApproximate - ApproxMatrix[:Q] * ApproxMatrix[:R]))

            #     # ApproxMatrix = psvdfact(OrigianlMatrixToApproximate, rank=Ranks[BlockIndex])
            #     # Error = maximum(abs.(OrigianlMatrixToApproximate -ApproxMatrix[:U] * Diagonal(ApproxMatrix[:S]) * ApproxMatrix[:Vt]))
                


            #     # println(Error)
            #     # println(Error2)
            #     # println(Ranks[BlockIndex])
            # end
            # # push!(ShearStiffness_H2,pqrfact(OrigianlMatrixToApproximate, rank=Ranks[i]))
            ApproxMatrix = pqrfact(OrigianlMatrixToApproximate, atol = AllowedError)
            push!(ShearStiffness_H2,ApproxMatrix)
            Ranks[BlockIndex] = size(ApproxMatrix[:Q],2)
            # push!(ShearStiffness_H2,psvdfact(OrigianlMatrixToApproximate, atol = 1e4))
            
        else 
            push!(ShearStiffness_H2,StiffnessMatrixShear[ElementRange_SR[i,1]:ElementRange_SR[i,2],ElementRange_SR[i,3]:ElementRange_SR[i,4]])
        end

    end

    ShearStiffness_H2 = ShearStiffness_H2[2:end]


    # ApproxMatrix_M = ApproxMatrix[:Q] * ApproxMatrix[:R]
    # ShearStiffness_H2 = ShearStiffness_H2[2:end]
    # figure(2)
    # plot(OrigianlMatrixToApproximate)
    # plot(ApproxMatrix[:Q] * ApproxMatrix[:R] .- OrigianlMatrixToApproximate)


    # figure(3)
    # clf()
    # i= 2000
    # # ApproxMatrix = pqrfact(OrigianlMatrixToApproximate, rank=10)
    # # ApproxMatrix_M = ApproxMatrix[:Q] * ApproxMatrix[:R]
    # ApproxMatrix = psvdfact(OrigianlMatrixToApproximate, rank=15)
    # ApproxMatrix_M = ApproxMatrix[:U] * Diagonal(ApproxMatrix[:S]) * ApproxMatrix[:Vt]
    # println(maximum(abs.(OrigianlMatrixToApproximate-ApproxMatrix_M)))
    # plot(ApproxMatrix_M[i,:])
    # plot(OrigianlMatrixToApproximate[i,:])

    # # 
    # # plot(ApproxMatrix_M[i,:] - OrigianlMatrixToApproximate[i,:])
    
    
    
    
    ThreadCount = 24
    BlockCount = length(Ranks)
    Par_ElementDivision2 = ParallelOptimization(ShearStiffness_H2, ElementRange_SR, FaultCount, BlockCount, ThreadCount)
    

    @time Elastic_Load_DispP= 
     HmatSolver_Pararllel(NetDisp, ShearStiffness_H2, ElementRange_SR, FaultCount, Par_ElementDivision2, ThreadCount)  
 
    @time ElasticLoadDisp = StiffnessMatrixShear * NetDisp
       
    maximum(abs.((ElasticLoadDisp - Elastic_Load_DispP)))

    @time for i=1:10
         Elastic_Load_DispP= 
     HmatSolver_Pararllel(NetDisp, ShearStiffness_H2, ElementRange_SR, FaultCount, Par_ElementDivision2, ThreadCount)  
    end

    @time for i=1:10
         ElasticLoadDisp = StiffnessMatrixShear * NetDisp
    end


    plot(ElasticLoadDisp)
    plot(Elastic_Load_DispP)
    ShearStiffness_H[1][:Q] * ShearStiffness_H[1][:R] - StiffnessMatrixShear[ElementRange_SR[1,1]:ElementRange_SR[1,2],ElementRange_SR[1,3]:ElementRange_SR[1,4]]
i=1
maximum(StiffnessMatrixShear[ElementRange_SR[i,1]:ElementRange_SR[i,2],ElementRange_SR[i,3]:ElementRange_SR[i,4]] - ShearStiffness_H2[1][:Q] * ShearStiffness_H2[1][:R])
ShearStiffness_H2[1] - StiffnessMatrixShear[ElementRange_SR[i,1]:ElementRange_SR[i,2],ElementRange_SR[i,3]:ElementRange_SR[i,4]]
# maximum(StiffnessMatrixShear[ElementRange_SR[i,1]:ElementRange_SR[i,2],ElementRange_SR[i,3]:ElementRange_SR[i,4]] - ShearStiffness_H2[1][:U] * ShearStiffness_H2[1][:Vt])

# @time for i=1:10
#     # @time Elastic_Load_Disp  = HmatSolver(NetDisp, ShearStiffness_H, BlockCount, ElementRange_SR, FaultCount)
#     @time Elastic_Load_DispP= 
#         HmatSolver_Pararllel(NetDisp, ShearStiffness_H, ElementRange_SR, FaultCount, Par_ElementDivision, ThreadCount) ./ Elastic_Load_DispP
#     @time ElasticLoadDisp = StiffnessMatrixShear * NetDisp
#         # println("One Step Over")
# end


# Elastic_Load_DispP

# println(Par_ElementDivision)
# SizeAndRankCumNorm = zeros(size(ShearStiffness_H,1))
# SizeAndRank = zeros(size(ShearStiffness_H,1))
# for i in eachindex(ShearStiffness_H[:])
#     if Ranks[i]>0
#         SizeAndRank[i] = size(ShearStiffness_H[i],1) * Ranks[i]
#     else 
#         SizeAndRank[i] = size(ShearStiffness_H[i],1)^2 
#     end

#     if i==1; SizeAndRankCumNorm[i] = SizeAndRank[i]
#     else
#         SizeAndRankCumNorm[i] =  SizeAndRankCumNorm[i-1] + SizeAndRank[i]
#     end
# end

# SizeAndRankCumNorm= SizeAndRankCumNorm/maximum(SizeAndRankCumNorm)


# # RecordedTime = ParallelOptimization(ShearStiffness_H, ElementRange_SR, FaultCount, BlockCount)
# CumRecordTime = cumsum(RecordedTime)
# CumRecordTimeNorm = CumRecordTime / maximum(CumRecordTime)
# # plot(CumRecordTimeNorm)

# ThreadCount = 24
# ElementGap = 1 / ThreadCount
# ElementCut = zeros(Int, ThreadCount +1)
# for i=1:ThreadCount
#     if i==ThreadCount 
#         ElementCut[i+1] = BlockCount
#     else 
#         ElementCut[i+1] = findmin(abs.(CumRecordTimeNorm .- ElementGap * i))[2]
#     end
# end
# Timeone = 0
# TimeP = 0

# println(ElementCut)
# println(Timeone)
# println(TimeP)
    
# findmin(abs.(SizeAndRankCumNorm .- 0.1))[2]


# plot(SizeAndRankCumNorm)
# ElementCut = zeros(Int,ThreadCount+1)
# for i=1:ThreadCount
#     if i==ThreadCount 
#         ElementCut[i+1] = BlockCount
#     else 
#         ElementCut[i+1] = ElementCut[i]  + round(10 * 1.4^i)
#     end
# end




# ElementCut = [0,20,60,90,350,1250,BlockCount]

# ElementCut = [0,10,20,35,50,65,85,900,BlockCount]
# ElementCut = [0,10,20,35,50,65,85,900,BlockCount]
# ElementCut = [0,2,40,60,100,350,550,1550,2800,BlockCount]
# ElementCut = [0,20,50,90,150,500,1000,2000,BlockCount]

# ThreadCount = length(ElementCut) - 1 

# println(maximum(abs.(Elastic_Load_Disp - Elastic_Load_DispP)))

# findmin(abs.(ElementCut .- 250))[2]
# Elastic_Load_Disp - Elastic_Load_DispP


# plot(Elastic_Load_Disp - Elastic_Load_DispP)


# Threads.nthreads()
# for i=1:10
#    println(size(ShearStiffness_H[i],1))
# end

# SizeDistribution = zeros(size(ShearStiffness_H,1))
# for i=1:size(ShearStiffness_H,1)
#     if i==1; SizeDistribution[i] = size(ShearStiffness_H[i],1);
#     else
#         SizeDistribution[i] =  SizeDistribution[i-1] + size(ShearStiffness_H[i],1)
#     end
# # size(ShearStiffness_H[1],1)

# end
# plot(SizeDistribution)
# Elastic_Load_Disp - Elastic_Load_DispP

#  for i=1:2
#   
# end

# @time Elastic_Load_Disp{2, Elastic_Load_Disp2, Elastic_Load_Disp3,Elastic_Load_Disp4} = 
#     HmatSolver_Pararllel(NetDisp, ShearStiffness_H, BlockCount, ElementRange_SR, FaultCount)
#      Elastic_Load_DispP = Elastic_Load_Disp1 + Elastic_Load_Disp2 + Elastic_Load_Disp3 +Elastic_Load_Disp4

# Elastic_Load_DispP - Elastic_Load_Disp
# ShearStiffness_H, ElasticLoadingShearMatrix, ReducedStiffnessMatrixShear, ReducedStiffnessMatrixNormal, Ranks, ElementRange_SR = 
#     load("HmatSave.jld2", "ShearStiffness_H", "ElasticLoadingShearMatrix", "ReducedStiffnessMatrixShear", "ReducedStiffnessMatrixNormal",
#          "Ranks", "ElementRange_SR")
# # FaultCount = length(StiffnessMatrixShear[:,1])

# NetDisp = ones(FaultCount)
# Elastic_Load_Disp = zeros(FaultCount)
# for Blockidx in 1:BlockCount             
#     Elastic_Load_Disp[ElementRange_SR[Blockidx,1]:ElementRange_SR[Blockidx,2]] = 
#         Elastic_Load_Disp[ElementRange_SR[Blockidx,1]:ElementRange_SR[Blockidx,2]] + 
#         ShearStiffness_H[Blockidx] * NetDisp[ElementRange_SR[Blockidx,3]:ElementRange_SR[Blockidx,4]]
# end

# @time Elastic_Load_DispOriginal=ElasticLoadingShearMatrix * (NetDisp) 

# Elastic_Load_DispOriginal-Elastic_Load_Disp

# # Ranks = RanksOriginal

=#