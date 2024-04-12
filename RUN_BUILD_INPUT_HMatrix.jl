
using LowRankApprox
using DelimitedFiles
using Base
using PyPlot
using PyCall
using JLD2
using Clustering
using LinearAlgebra
@pyimport matplotlib.patches as patches
pygui(true)


include("Functions_BuildInputFile.jl")
include("Functions_OKADA3D.jl")
include("Results/Functions_Plot.jl")
include("Functions_Hmatrix.jl")



InputBulkFileName="Input_BulkFaultGeometry.txt"
OutputFileName="Input_Discretized.jld2"

ArrangePoint = [0,-3000,1000]
HowManyDivisionEachLevel = 2
TotalHierarchyLevel = 8
MinimumElementsToCut = 100
ElementPartRoughCount = 200
DistDiamRatioCrit = 1

############# Plots? ##############
PlotHMat = 1
PlotBlock3D = 1
###################################





############################# Read Bulk Input ##################################
######++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++######
#### Bulk File Order 
####  1.SwitchSSRN  2.ShearMod  3.PoiRat  4.R_Density   
####  5. Crit_TooClose   6. TooCloseNormal_Multiplier
####  ----------------------------------------------------------------------------
####  1.Ctr_X     2.Ctr_Y 3.Ctr_Z 4.St_L	    5.Dip_L	    6.StAng	    7.DipAng	8.LR/RN
####  9.a         10.b	11.Dc	12.Theta_i	13. V_i     14. Friction_i 15.NormalStress at surface [Pa]  
####  16. NoarmalStress Gradient [Pa] 17. V_Const     18. Minimum Segment Length

Input_Bulk=readdlm(InputBulkFileName)

Switch_StrikeSlip_or_ReverseNormal = Input_Bulk[2,1] 
ShearModulus = Input_Bulk[2,2]
PoissonRatio = Input_Bulk[2,3]
RockDensity = Input_Bulk[2,4]
DropCrit= Input_Bulk[2,5]
DropCritNormalStressMultiplier= Input_Bulk[2,6]
MinimumNS=Input_Bulk[2,7]

Input_Bulk=Input_Bulk[4:end,:]
Input_Bulk=Input_Bulk[sortperm(Input_Bulk[:, 17]), :]


# Adjust if positive depth exists
# for i=1:size(Input_Bulk,1)
for i in eachindex(Input_Bulk[:,1])
    if Input_Bulk[i,3] < Input_Bulk[i,5] / 2*sind(Input_Bulk[i,7])
        println("Caution! Fault ",i," may have negative depth")
    end
end

########^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^########
################################################################################






############################### Bulk to Segment ################################
######++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++######
#### Segment File Order 
####  1.Ctr_X     2.Ctr_Y 3.Ctr_Z 4.St_L	    5.Dip_L	    6.StAng	    7.DipAng	8.LR/RN
####  9.a         10.b	11.Dc	12.Theta_i	13. V_i     14. Friction_i 15.NormalStress  
####  16. V_Const 17. Bulk Index     18. Bulk Strike Length      19. Bulk Dip Length


Input_Segment = BulkToSegment(Input_Bulk);
FaultCount=   size(Input_Segment,1)
Input_Segment=[Input_Segment ones(FaultCount) zeros(FaultCount)]



################################# Group and Sort #################################
Block_Ctr_Diam, Block_Range_Level, Input_Segment = 
    GroupAndSort_AllLevel(HowManyDivisionEachLevel, TotalHierarchyLevel, MinimumElementsToCut,
                         ArrangePoint, FaultCount)
println("Grouping and Sorting Done")


############################## Build Hierarchy ####################################
Ranks, ElementRange_SR = BuildHierarchy(Block_Range_Level, Block_Ctr_Diam, DistDiamRatioCrit) 




################################################################################
#################################### 3D Plot ###################################

if PlotBlock3D == 1
    println("Plotting 3D Group")
    PlotRotation=[45,-30]
    Edge = 0
    Transparent = 0
    MinMax_Axis=0 # automatically detect max and min 
    LoadingFaultCount=0 
    ColorMinMax = 0  
    figure(1)
    clf()
    MaxVaule, MinValue = FaultPlot_3D_Color_General(Input_Segment[:,1:3],Input_Segment[:,4], Input_Segment[:,5],
        Input_Segment[:,6], Input_Segment[:,7], Input_Segment[:,8], Input_Segment[:,20], 
        PlotRotation, MinMax_Axis, ColorMinMax, Transparent, Edge, LoadingFaultCount)

    # figure(1)
    # plotforcbar=  scatter([1,1],[1,1],0.1, [MinValue,MaxVaule], cmap="jet")
    plotforcbar=  scatter([1,1],[1,1],0.1, [MinValue,MaxVaule], cmap="jet")
    colorbar(plotforcbar, pad=0.15)
    figure(1).canvas.draw()
    xlabel("x")
    ylabel("y")
end
########^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^########
################################################################################



################################################################################
########################### Plot Hierarchical Matrix ###########################
if PlotHMat == 1
    figure(2)
    println("Plotting Hmatrix Structure")
    clf()
    ax = gca()
    # ax[:set_aspect]("equal")
    for i=1:length(ElementRange_SR[:,1])
        if Ranks[i] > 6 
            c = PyObject(patches.Rectangle((ElementRange_SR[i,3], -ElementRange_SR[i,1]), ElementRange_SR[i,4] - ElementRange_SR[i,3], 
                        -ElementRange_SR[i,2] + ElementRange_SR[i,1], linewidth=1, edgecolor="k", facecolor=[0.2 0.3 1]))                                
        elseif Ranks[i] > 0                 
            c = PyObject(patches.Rectangle((ElementRange_SR[i,3], -ElementRange_SR[i,1]), ElementRange_SR[i,4] - ElementRange_SR[i,3], 
                        -ElementRange_SR[i,2] + ElementRange_SR[i,1], linewidth=1, edgecolor="k", facecolor=[0.5 0.8 1]))
        else           
            c = PyObject(patches.Rectangle((ElementRange_SR[i,3], -ElementRange_SR[i,1]), ElementRange_SR[i,4] - ElementRange_SR[i,3], 
                        -ElementRange_SR[i,2] + ElementRange_SR[i,1], linewidth=1, edgecolor="k", facecolor=[1 0.8 0.4]))
        end
        ax.text( (ElementRange_SR[i,3] + ElementRange_SR[i,4])/2, -(ElementRange_SR[i,1] + ElementRange_SR[i,2])/2, i ,size=8, horizontalalignment="center", verticalalignment="center", color="k") 
        ax.add_patch(c) 
    end
    xlim(0,FaultCount)
    ylim(-FaultCount,0)
end
########^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^########
################################################################################


    
################### Build Matrix by Parts #################
StiffnessMatrixShearOriginal, StiffnessMatrixNormalOriginal = 
    BuildMatrixByPartsShear(FaultCount, ElementPartRoughCount, Input_Segment,  ShearModulus, PoissonRatio)

jldsave("HmatSave.jld2"; StiffnessMatrixShearOriginal, StiffnessMatrixNormalOriginal, Ranks, ElementRange_SR)

#=

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
        # push!(ABlock_H,psvdfact(ShearMatrix[ElementRange_SR[i,1]:ElementRange_SR[i,2],ElementRange_SR[i,3]:ElementRange_SR[i,4]], rank=Rank[i]))
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
# for i=1:10
@time begin
    for i=1:BlockCount
        AX_Approx[ElementRange_SR[i,1]:ElementRange_SR[i,2]] = 
           AX_Approx[ElementRange_SR[i,1]:ElementRange_SR[i,2]] + ABlock_H[i] * X[ElementRange_SR[i,3]:ElementRange_SR[i,4]]
    end
end

@time AX = ShearMatrix * X

AX_Error = AX - AX_Approx
println(maximum(abs.(AX_Error))/maximum(abs.(AX)))

figure(3)
clf()
plot(AX)
plot(AX_Approx)



# for i=1:BlockedElementCount
#     for j=1:BlockedElementCount
#         ElementIdx=ElementIdx+1
#         ElementRange_SR[ElementIdx,:] = 
#         [minimum(findall(x->x==i, Input_Segment[:,20])) , maximum(findall(x->x==i, Input_Segment[:,20])) ,
#         minimum(findall(x->x==j, Input_Segment[:,20])) , maximum(findall(x->x==j, Input_Segment[:,20]))]
#     end
# end



#=


############### Admissibility and Element Set ##################

ElementCount::Int= length(Input_Segment[:,20])
BlockedElementCount::Int = maximum(Input_Segment[:,20])
BlockCount = BlockedElementCount^2
Ranks =  zeros(Int,BlockCount)

BlockCenter = zeros(3, BlockedElementCount)
BlockDiameter = zeros(BlockedElementCount)
BlockDistance = zeros(BlockedElementCount,BlockedElementCount)
BlockCenter = ClusterCenter
BlockDiameter = 1000*ones(FaultCount)

# for i=1:BlockedElementCount
#     BlockCenter[:,i] = ClusterCenter[:,GroupOrder[i]]
#     BlockDiameter[i] = 1000
# end


BlockIndex = 0
for i=1:BlockedElementCount
    for j=1:BlockedElementCount
        BlockIndex = BlockIndex + 1 
        BlockDistance[i,j] = norm(BlockCenter[:,i] - BlockCenter[:,j])
        if (BlockDiameter[i] + BlockDiameter[j])/2 < BlockDistance[i,j] 
            # println(i,"  ",j)
            Ranks[BlockIndex] = 2
        end
    end
end


ElementRange_SR = zeros(Int, BlockCount,4)
ElementIdx=0
for i=1:BlockedElementCount
    for j=1:BlockedElementCount
        ElementIdx=ElementIdx+1
        ElementRange_SR[ElementIdx,:] = 
        [minimum(findall(x->x==i, Input_Segment[:,20])) , maximum(findall(x->x==i, Input_Segment[:,20])) ,
        minimum(findall(x->x==j, Input_Segment[:,20])) , maximum(findall(x->x==j, Input_Segment[:,20]))]
    end
end







################################################################################
#################################### 3D Plot ###################################
PlotRotation=[45,-30]
Edge = 0
Transparent = 0
MinMax_Axis=0 # automatically detect max and min 
LoadingFaultCount=0 
ColorMinMax = 0  
figure(1)
clf()
MaxVaule, MinValue = FaultPlot_3D_Color_General(Input_Segment[:,1:3],Input_Segment[:,4], Input_Segment[:,5],
    Input_Segment[:,6], Input_Segment[:,7], Input_Segment[:,8], Input_Segment[:,20], 
    PlotRotation, MinMax_Axis, ColorMinMax, Transparent, Edge, LoadingFaultCount)

# figure(1)
plotforcbar=  scatter([1,1],[1,1],0.1, [MinValue,MaxVaule], cmap="jet")
colorbar(plotforcbar, pad=0.15)
figure(1).canvas.draw()
xlabel("x")
ylabel("y")
########^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^########
################################################################################

    
# # # ##################### Build Matrix #####################

# StiffnessMatrixShearOriginal, StiffnessMatrixNormalOriginal =
# Function_StiffnessMatrixStrikeSlip_Vector(Input_Segment, ShearModulus, PoissonRatio)

# ShearMatrix = StiffnessMatrixShearOriginal



# ##################### Approximate #####################
# ABlock = Any[0]
# ABlock_H = Any[0]
# MatrixElement = 0
# for i=1:BlockCount
#     MatrixElement = MatrixElement + 1
#     println(MatrixElement)
    
#     # push!(ABlock_H,A[Block_Edges[i,1]+1:Block_Edges[i,2],Block_Edges[i,3]+1:Block_Edges[i,4]])

#     if Ranks[MatrixElement] > 0
#         # push!(ABlock_H,psvdfact(ShearMatrix[ElementRange_SR[i,1]:ElementRange_SR[i,2],ElementRange_SR[i,3]:ElementRange_SR[i,4]], rank=Rank[i]))
#         push!(ABlock_H,pqrfact(ShearMatrix[ElementRange_SR[i,1]:ElementRange_SR[i,2],ElementRange_SR[i,3]:ElementRange_SR[i,4]], rank=Ranks[i]))
#     else 
#         push!(ABlock_H,ShearMatrix[ElementRange_SR[i,1]:ElementRange_SR[i,2],ElementRange_SR[i,3]:ElementRange_SR[i,4]])
#     end

# end
# ABlock_H = ABlock_H[2:end]


# ####################### Calculate ########################

# X = ones(FaultCount)
# # X = rand(n)
# # for i=1:10
# @time begin
#     AX_Approx = zeros(FaultCount)
#     for i=1:BlockCount
#         AX_Approx[ElementRange_SR[i,1]:ElementRange_SR[i,2]] = AX_Approx[ElementRange_SR[i,1]:ElementRange_SR[i,2]] + ABlock_H[i] * X[ElementRange_SR[i,3]:ElementRange_SR[i,4]]
#     end
# end
# # end

# @time AX = ShearMatrix * X
# AX_Error = AX - AX_Approx
# println(maximum(abs.(AX_Error)))

# figure(2)
# clf()
# plot(AX)
# plot(AX_Approx)



# # Input_Segment=sortslices(Input_Segment, dims = 1, by = x -> c)


# # # Input_Segment[1:800,:]=sortslices(Input_Segment[1:800,:], dims = 1, by = x -> c)



# # MinMaxX= [minimum(Input_Segment[:,1]), maximum(Input_Segment[:,1])]
# # MinMaxY= [minimum(Input_Segment[:,2]), maximum(Input_Segment[:,2])]
# # MinMaxZ= [minimum(Input_Segment[:,3]), maximum(Input_Segment[:,3])]
# # figure(2)
# # plot(FaultCenter[:,1]/maximum(abs.(FaultCenter[:,1])))
# # plot(a)

=#
=#