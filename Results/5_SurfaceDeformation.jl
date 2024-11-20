
using PyPlot
using PyCall
using JLD2
using LinearAlgebra
using Statistics
pygui(true)
include("Functions_Plot.jl")
include("../Functions_OKADA3D.jl")

ResultName="Result"
FileName="Results/" * ResultName * ".jld2"
FileNameInput="Results/" * ResultName * "_Input.jld2"

ResultTime, ResultV, ResultDisp, Result_NormalStress, History_Theta, =
load(FileName,"History_Time", "History_V", "History_Disp", "History_NormalStress","History_Theta")
ResultV[ResultV.<=0] .= 1e-100

#######################################################################################
############################### Figure Configuration ##################################

# figure(10); clf(); PyPlot.plot(log10.(ResultV[:,1:30:end])); xlabel("Record Step")
PlotStep = 500
SurfacePlotRangeX = [-5000, 30000]
SurfacePlotRangeY = [-5000, 40000]
SurfacePlotIntervalX = 1000
SurfacePlotIntervalY = 1000

################## What to Plot ? ###################
# PlotInput=log10.(ResultV[PlotStep,:]); ColorMinMax=[-12,0] 
PlotInput=ResultDisp[PlotStep,:]  ; ColorMinMax=0 
#############---------------------------#############


############# Saving Multiple Figures ###############
Animation_Save = 0 # 1 for save
StepBegin = 10 # first record step
StepEnd = 100
StepInterval = 10
###################^^^^^^^^^^^^^^^###################

##############+++++++++++++++++++++++++++++++++++++++++++++++++++++++++################
#######################################################################################

FaultCount= load(FileNameInput, "FaultCount")
FaultCenter= load(FileNameInput, "FaultCenter")
FaultLengthStrike= load(FileNameInput, "FaultLengthStrike")
FaultLengthDip= load(FileNameInput, "FaultLengthDip")
FaultStrikeAngle= load(FileNameInput, "FaultStrikeAngle")
FaultDipAngle= load(FileNameInput, "FaultDipAngle")
FaultLLRR= load(FileNameInput, "FaultLLRR")
Fault_BulkIndex= load(FileNameInput, "Fault_BulkIndex")
ShearModulus= load(FileNameInput, "ShearModulus")
PoissonRatio= load(FileNameInput, "PoissonRatio")



XElementCount = Int(round((SurfacePlotRangeX[2] - SurfacePlotRangeX[1]) / SurfacePlotIntervalX +1))
SurfacePlotIntervalAdjX = (SurfacePlotRangeX[2] - SurfacePlotRangeX[1]) / (XElementCount-1)
YElementCount = Int(round((SurfacePlotRangeY[2] - SurfacePlotRangeY[1]) / SurfacePlotIntervalY +1))
SurfacePlotIntervalAdjY = (SurfacePlotRangeY[2] - SurfacePlotRangeY[1]) / (YElementCount-1)

PlotXCountour = zeros(XElementCount,YElementCount)
PlotYCountour = zeros(XElementCount,YElementCount)
PlotSurface = zeros(XElementCount,YElementCount)
DistanceSum = zeros(XElementCount,YElementCount)
ReceiverCenter = zeros(XElementCount * YElementCount, 3)

receiveridx = 0
for xidx = 1:XElementCount
    for yidx = 1:YElementCount
    receiveridx += 1
        XYZ = [SurfacePlotRangeX[1] + SurfacePlotIntervalAdjX * (xidx-1), SurfacePlotRangeY[1] + SurfacePlotIntervalAdjY * (yidx-1), 0]
        PlotXCountour[xidx,yidx] = XYZ[1]
        PlotYCountour[xidx,yidx] = XYZ[2]
        ReceiverCenter[receiveridx,:] = [XYZ[1], XYZ[2], 0]
    end
end



# ReceiverCenter = [0 0 0]


FaultCountSource = FaultCount
FaultCenterSource = FaultCenter
FaultLengthStrikeSource = FaultLengthStrike
FaultLengthDipSource =  FaultLengthDip
FaultStrikeAngleSource = FaultStrikeAngle
FaultDipAngleSource = FaultDipAngle
FaultLLRRSource = FaultLLRR  
FaultCenterReceiver = ReceiverCenter


FaultCountReceiver=length(FaultCenterReceiver[:,1])
# FaultCenterReceiver=Input_SegmentReceiver[:,1:3]
FaultStrikeAngleReceiver= zeros(FaultCountReceiver)
FaultDipAngleReceiver = zeros(FaultCountReceiver)
FaultLLRRReceiver = ones(FaultCountReceiver)
# println(FaultCountSource, "  ", FaultCountReceiver)

DisplacementX = zeros(FaultCountReceiver)
DisplacementY = zeros(FaultCountReceiver)
DisplacementZ = zeros(FaultCountReceiver)

for SourceIndex=1:FaultCount

    # SourceIndex=20

    println(SourceIndex,"  ",SourceIndex,"/",FaultCountSource)

    SourceCenter = FaultCenterSource[SourceIndex,:];
    SourceLengthStrike = FaultLengthStrikeSource[SourceIndex];
    SourceLengthDip = FaultLengthDipSource[SourceIndex];
    SourceStrikeAngle = FaultStrikeAngleSource[SourceIndex];
    SourceDipAngle = FaultDipAngleSource[SourceIndex];
    SourceLLRR = FaultLLRRSource[SourceIndex];
            
    ReceiverStrikeAngle = FaultStrikeAngleReceiver
    ReceiverDipAngle = FaultDipAngleReceiver
    ReceiverLLRR = FaultLLRRReceiver
    RelativeStrkieAngle = ReceiverStrikeAngle .- SourceStrikeAngle

    # DISL1 = -SourceLLRR; # Left Latteral is +1 for Okada
       
    DISL1 = -PlotInput[SourceIndex]
    DISL2 = 0.0;
    DISL3 = 0.0;            
            
    DEPTH=SourceCenter[3] # Source Depth
    AL1=SourceLengthStrike/2
    AL2=SourceLengthStrike/2
    AW1=SourceLengthDip/2
    AW2=SourceLengthDip/2

    LameFirstParam=2*ShearModulus*PoissonRatio/(1-2*PoissonRatio);
    ALPHA=(LameFirstParam+ShearModulus)/(LameFirstParam+2*ShearModulus);

    #######################################################################
    ##### Calculate Receiver Point Relative to the Source and Source frame

    X_Dist = ReceiverCenter[:,1] .- SourceCenter[1]
    Y_Dist = ReceiverCenter[:,2] .- SourceCenter[2]

    X = X_Dist .* cosd(-SourceStrikeAngle) .- Y_Dist .* sind(-SourceStrikeAngle)
    Y = X_Dist .* sind(-SourceStrikeAngle) .+ Y_Dist .* cosd(-SourceStrikeAngle)
    Z = zeros(FaultCountReceiver)

    #######################################################################
    ##### Calculate Stress Change at Source Frame

    UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET = Okada_DC3D_Vector(ALPHA,
        X,Y,Z,DEPTH,SourceDipAngle,
        AL1,AL2,AW1,AW2,DISL1,DISL2,DISL3);


    for ReceiverIdx = 1:FaultCountReceiver
        
        ReceiverDipAngle = 0.0
        ReceiverLLRR = 1.0
        RelativeStrkieAngle = - SourceStrikeAngle

        UXYZ_SourceFrame=[UX[ReceiverIdx], UY[ReceiverIdx], UZ[ReceiverIdx]]


        #######################################################################
        ##### Rotate Source Frame Stress to Flat Receiver 

        RotationMat_FromReceiver_Strike=
        [cosd(-RelativeStrkieAngle) -sind(-RelativeStrkieAngle)  0
        sind(-RelativeStrkieAngle) cosd(-RelativeStrkieAngle) 0
        0  0  1];

                                
        UXYZ_Receiver = RotationMat_FromReceiver_Strike * UXYZ_SourceFrame

        UXYZ_Receiver = RotationMat_FromReceiver_Strike * UXYZ_SourceFrame
        if sum(isnan.(UXYZ_Receiver))>0; println("nan here"); end
        UXYZ_Receiver[isnan.(UXYZ_Receiver)] .= 0.0

        DisplacementX[ReceiverIdx] += UXYZ_Receiver[1]
        DisplacementY[ReceiverIdx] += UXYZ_Receiver[2]
        DisplacementZ[ReceiverIdx] += UXYZ_Receiver[3]

    end

end






receiveridx = 0
for xidx = 1:XElementCount
    for yidx = 1:YElementCount
        receiveridx += 1
        PlotSurface[xidx,yidx] = DisplacementY[receiveridx]
        # ReceiverCenter[receiveridx,:] = [XYZ[1], XYZ[2], 0]
    end
end






figure(1)
clf()
contourf(PlotXCountour,PlotYCountour,-PlotSurface, cmap= "jet", 100,vmin=-0.3, vmax=0.3)
# contourf(PlotXCountour,PlotYCountour,PlotSurface, cmap= "jet", 100,line_smoothing=5)
# plot_surface(PlotXCountour,PlotYCountour,PlotSurface, 20)

colorbar()
# ContourX = zeros
