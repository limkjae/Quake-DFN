
using DelimitedFiles
using PyPlot
using PyCall

pygui(true)
include("Functions_BuildInputFile.jl")
include("Results/Functions_Plot.jl")
InputBulkFileName="Input_BulkFaultGeometry.txt"
Input_Bulk=readdlm(InputBulkFileName)
RakeLRRN = Int(Input_Bulk[2,1])
Input_Bulk=Input_Bulk[4:end,:]
Faultcount = length(Input_Bulk[:,1])


#######################################################################################
############################### Figure Configuration ##################################

PlotRotation=[30,-50]
Transparent = 1 # 1 for transparent fault plot
Edge = 1 # 0 for no element boudary 
MinMax_Axis=0
LoadingFaultCount=0 # 1 to plot constant velocity faults. 

#######################################################################################



RakeAngle = Input_Bulk[:,8]
DipAngle =  Input_Bulk[:,7]
# adjust sense of slip to rake angle
if RakeLRRN ==1
    for BulkIndex = 1: Faultcount
        if RakeAngle[BulkIndex] == -1.0
            RakeAngle[BulkIndex] = 0.0
        else        
            RakeAngle[BulkIndex] = 180.0
        end    
    end
elseif RakeLRRN ==2
    for BulkIndex = 1: Faultcount
        if DipAngle[BulkIndex] < 90.0
            if RakeAngle[BulkIndex] == -1.0
                RakeAngle[BulkIndex] = 90.0
            else        
                RakeAngle[BulkIndex] = 270.0
            end    
        else 
            if RakeAngle[BulkIndex] == -1.0
                RakeAngle[BulkIndex] = 270.0
            else        
                RakeAngle[BulkIndex] = 90.0
            end    
        end
    end
end

RakeAngle_NRAdjusted = copy(RakeAngle)
for BulkIndex = 1: Faultcount
    if DipAngle[BulkIndex] > 90.0
        RakeAngle_NRAdjusted[BulkIndex] = 360 - RakeAngle_NRAdjusted[BulkIndex]

    end
end
PlotInput=RakeAngle_NRAdjusted; ColorMinMax=[0, 360]

    


figure(8)
clf()
ColorMap = "jet"
MaxVaule, MinValue = FaultPlot_3D_Color_General_hsv(Input_Bulk[:,1:3],
    Input_Bulk[:,4], Input_Bulk[:,5], Input_Bulk[:,6], Input_Bulk[:,7], Input_Bulk[:,8], PlotInput, 
    PlotRotation, MinMax_Axis, ColorMinMax, Transparent, Edge, LoadingFaultCount)
    ax = subplot(projection="3d")
    xlabel("x")
    ylabel("y")
plotforcbar=  scatter([1,1],[1,1],0.1, [MinValue,MaxVaule], cmap="hsv")
cbar  = colorbar(plotforcbar, pad=0.15)
cbar.set_ticks([0, 90, 180,270,360])
cbar.set_ticklabels(["Left Latteral", "Reverse", "Right Lateral", "Normal", "Left Lateral"])
figure(8).canvas.draw()






# PlotInput=log10.(Input_Bulk[:,13]); ColorMinMax=0

##### Input Code
##### 1.Ctr_X     2.Ctr_Y 3.Ctr_Z 4.St_L	    5.Dip_L	    6.StAng	    7.DipAng	8.LR/RN
#####  9.a         10.b	    11.Dc	12.Theta_i	13. V_i     14. Friction_i 15.NormalStress at surface [Pa]  
#####  16. NoarmalStress Gradient [Pa] 17. V_Const     18. Minimum Segment Length


LineLength = vec(minimum([Input_Bulk[:,4] Input_Bulk[:,5]], dims=2)./2)

UnrotatedSlipUnitVec = [cosd.(RakeAngle) sind.(RakeAngle) zeros(Faultcount)]
UnrotatedGapVector = [0 0 1]
RotatedSlipUnitVec = UnrotatedSlipUnitVec .* 0.0 
VecStart = UnrotatedSlipUnitVec .* 0.0 
VecEnd = UnrotatedSlipUnitVec .* 0.0 
FaultCenter =  Input_Bulk[:,1:3]
FaultCenter[:,3] = -FaultCenter[:,3]
for BulkIndex = 1: Faultcount

RotationMat_Strike=
[cosd(Input_Bulk[BulkIndex,6]) -sind(Input_Bulk[BulkIndex,6])  0
sind(Input_Bulk[BulkIndex,6]) cosd(Input_Bulk[BulkIndex,6]) 0
0  0  1];

RotationMat_Dip=
[1 0 0
0 cosd(Input_Bulk[BulkIndex,7]) -sind(Input_Bulk[BulkIndex,7])
0 sind(Input_Bulk[BulkIndex,7]) cosd(Input_Bulk[BulkIndex,7])]

RotatedSlipUnitVec[BulkIndex,:] = RotationMat_Strike * RotationMat_Dip  * UnrotatedSlipUnitVec[BulkIndex,:]
RotatedGap = RotationMat_Strike * RotationMat_Dip  * UnrotatedGapVector' .* 50

VecStart[BulkIndex,:] = -RotatedSlipUnitVec[BulkIndex,:] /2 * LineLength[BulkIndex] .+ FaultCenter[BulkIndex,1:3] + RotatedGap
VecEnd[BulkIndex,:] = RotatedSlipUnitVec[BulkIndex,:]/2* LineLength[BulkIndex] .+ FaultCenter[BulkIndex,1:3] - + RotatedGap
end



# fig = plt.figure()
# ax = subplot(projection="3d")
for i =1:Faultcount 
    # ax.plot([VecStart[i,1], VecEnd[i,1]], [VecStart[i,2],VecEnd[i,2]],zs=[VecStart[i,3],VecEnd[i,3]],"k")
    ax.quiver(VecStart[i,1], VecStart[i,2], VecStart[i,3], 
        RotatedSlipUnitVec[i,1] * LineLength[i], RotatedSlipUnitVec[i,2] * LineLength[i], RotatedSlipUnitVec[i,3] * LineLength[i],
        color="k",arrow_length_ratio=0.2)
    ax.quiver(VecEnd[i,1], VecEnd[i,2], VecEnd[i,3], 
        -RotatedSlipUnitVec[i,1] * LineLength[i], -RotatedSlipUnitVec[i,2] * LineLength[i], -RotatedSlipUnitVec[i,3] * LineLength[i],
        color="k",arrow_length_ratio=0.2)

end



