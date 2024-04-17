using JLD2
using LowRankApprox
using LinearAlgebra
using PyPlot
using PyCall
using Distributed
pygui(true)

# jldsave("HmatSave.jld2"; StiffnessMatrixShearOriginal, StiffnessMatrixNormalOriginal, Ranks, ElementRange_SR)


ShearMatrix, NormalMatrix, Ranks, ElementRange_SR = load("HmatSave.jld2", "ReducedStiffnessMatrixShear", "ReducedStiffnessMatrixNormal", 
            "Ranks", "ElementRange_SR", "ShearStiffness_H", "ElasticLoadingShearMatrix")
FaultCount = length(ShearMatrix[:,1])

# Ranks = RanksOriginal


##################### Approximate #####################
BlockCount = length(ElementRange_SR[:,1])
ABlock = Any[0]
ABlock_H = Any[0]
BlockIndex = 0
# RanksOriginal = Ranks
Ranks = Ranks .* 2
for i=1:BlockCount
    BlockIndex = BlockIndex + 1
    # println(MatrixElement)
    
    # push!(ABlock_H,A[Block_Edges[i,1]+1:Block_Edges[i,2],Block_Edges[i,3]+1:Block_Edges[i,4]])

    if Ranks[BlockIndex] > 0
        # push!(ABlock_H,psvdfact(ShearMatrix[ElementRange_SR[i,1]:ElementRange_SR[i,2],ElementRange_SR[i,3]:ElementRange_SR[i,4]], rank=Ranks[i]))
        push!(ABlock_H,pqrfact(ShearMatrix[ElementRange_SR[i,1]:ElementRange_SR[i,2],ElementRange_SR[i,3]:ElementRange_SR[i,4]], rank=Ranks[i]))
    else 
        push!(ABlock_H,ShearMatrix[ElementRange_SR[i,1]:ElementRange_SR[i,2],ElementRange_SR[i,3]:ElementRange_SR[i,4]])
    end

end
ABlock_H = ABlock_H[2:end]

################################ Test ##############################
X = ones(FaultCount)
# X = rand(n)
# for i=1:10

AX_Approx = zeros(FaultCount)
# @time begin
    @distributed for i in 1:BlockCount        
    for i in 1:BlockCount
                AX_Approx[ElementRange_SR[i,1]:ElementRange_SR[i,2]] = 
                    AX_Approx[ElementRange_SR[i,1]:ElementRange_SR[i,2]] + ABlock_H[i] * X[ElementRange_SR[i,3]:ElementRange_SR[i,4]]
                
    end
end

@time AX = ShearMatrix * X


# end

AX_Error = AX - AX_Approx
println(maximum(abs.(AX_Error))/maximum(abs.(AX)))

figure(3)
clf()
plot(AX)
plot(AX_Approx)
