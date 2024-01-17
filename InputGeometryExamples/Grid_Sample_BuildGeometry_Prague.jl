using DelimitedFiles
using PyPlot
using PyCall
using LinearAlgebra
pygui(true)

include("../Functions_BuildInputFile.jl")


function BuildBulk()

    OutputFileName="Input_BulkFaultGeometry.txt"

    SwitchSSRN = 1
    ShearMod= 20e9
    PoissonRatio = 0.2
    R_Density=2400.0
    Crit_TooClose= 1.05
    TooCloseNormal_Multiplier = 0.6
    MinNormalStress = 2e6
    
    FaultDepthCenter=2000.0
    Multiple=110

    MaxSigSurface = 2e6
    MaxSigGrad = 10000.0 #Pa/m
    MaxStressOrientation = 10. # between 0-180 degree
    StressRatioMaxOverMin = 0.5
    MinFriction = 0.2

    RawData=readdlm("InputGeometryExamples/CustomFaultInfoOK.txt", '\t')

    X1=RawData[:,1]
    X2=RawData[:,1] + RawData[:,3]
    Y2=-RawData[:,2]
    Y1=-RawData[:,2] - RawData[:,4]
    for i=17:length(RawData[:,1])
        
        Y2Temp=Y2[i]
        Y2[i]=Y1[i]
        Y1[i]=Y2Temp
    end
    InjectorCenter=[94.74, 125.734]

    X1=(X1 .- InjectorCenter[1]) * Multiple
    Y2=(Y2 .+ InjectorCenter[2]) * Multiple
    X2=(X2 .- InjectorCenter[1]) * Multiple
    Y1=(Y1 .+ InjectorCenter[2]) * Multiple

    BulkFaultcount=length(X1)
    FaultDistacne = [X1 Y1] - [X2 Y2]
    FaultLength = zeros(BulkFaultcount)
    FaultStrikeAngle = zeros(BulkFaultcount)
    for i=1:BulkFaultcount
        FaultLength[i] = norm([FaultDistacne[i,1],FaultDistacne[i,2]])
        FaultStrikeAngle[i] = rad2deg(atan(FaultDistacne[i,2]/FaultDistacne[i,1]))
    end


    F_a = 0.003 
    F_b =0.006 
    F_Dc = 1e-3
    F_Theta = 1e8
    F_Vini = 1e-15
    MaximumSegLength=200 #200

    ############################# Write Bulk Input #################################
    ######++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++######
    ### Bulk File Order 
    ###  1.Ctr_X     2.Ctr_Y 3.Ctr_Z 4.St_L	    5.Dip_L	    6.StAng	    7.DipAng	8.LR/RN (-1: LL / +1 RL or -1:Reverse / +1 Normal)
    ###  9.a         10.b	11.Dc	12.Theta_i	13. V_i     14. Friction_i 15.NormalStress at surface [Pa]  
    ###  16. NoarmalStress Gradient [Pa] 17. V_Const     18. Minimum Segment Length


    FaultInfo=zeros(BulkFaultcount,18)
    FaultInfo[:,1] = (X1 + X2)/2
    FaultInfo[:,2] = (Y1 + Y2)/2
    FaultInfo[:,3] .= FaultDepthCenter;
    FaultInfo[:,4] = FaultLength
    FaultInfo[:,5] .= FaultDepthCenter*2
    FaultInfo[:,6] = FaultStrikeAngle;
    FaultInfo[:,7] .= 90.0;
    FaultInfo[:,9] .= F_a;
    FaultInfo[:,10] .= F_b;
    FaultInfo[:,11] .= F_Dc;
    FaultInfo[:,12] .= F_Theta;
    FaultInfo[:,13] .= F_Vini;
    FaultInfo[:,17] .= 0.0;
    FaultInfo[:,18] .= MaximumSegLength;


    ### Friction and NormalSTress
    NormalStressParameter = (1+StressRatioMaxOverMin)/2 .+ (1-StressRatioMaxOverMin)/2 .* cosd.(2 .* (FaultInfo[:,6] .- 90.0 .- MaxStressOrientation))
    ShearStressParameter = (1-StressRatioMaxOverMin)/2 .* sind.(2 .* (FaultInfo[:,6] .- 90.0 .- MaxStressOrientation))
    FaultInfo[:,8].= sign.(ShearStressParameter);
    FaultInfo[:,14] = abs.(ShearStressParameter ./ NormalStressParameter);
    FaultInfo[:,15] =MaxSigSurface.* NormalStressParameter;
    FaultInfo[:,16] =MaxSigGrad .* NormalStressParameter;
    
    for i=1:length(FaultInfo[:,8])
        if FaultInfo[i,8] == 0
            FaultInfo[i,8] = 1
        end
    end
    
    for i=1:length(FaultInfo[:,14])
        if FaultInfo[i,14] < MinFriction
            FaultInfo[i,14] = MinFriction
        end
    end

    figure(1)
    clf()
    FaultPlot_3D(FaultInfo[:,1:3],FaultInfo[:,4], FaultInfo[:,5], FaultInfo[:,6], FaultInfo[:,7], FaultInfo[:,8])
    xlabel("x")
    ylabel("y")



    InputFile=[FaultInfo;]

    FaultSegmentCount=zeros(1,2)
    for i in eachindex(InputFile[:,1])    
        FaultSegmentCount=[FaultSegmentCount;[ceil(InputFile[i,4]/InputFile[i,18]), ceil(InputFile[i,5]/InputFile[i,18])]']
    end
    FaultSegmentCount=FaultSegmentCount[2:end,:]
    TotalSegmentCount=sum(Int64,FaultSegmentCount[:,1].*FaultSegmentCount[:,2])
    println("TotalSegment Count will be ", TotalSegmentCount)


    open(OutputFileName, "w") do io
        write(io,"SwitchSS/RN\tShearMod\tPoissonRatio\tR_Density\tCrit_TooClose\tTooCloseNormal_Multiplier\tMinimum_NS\n")
        writedlm(io,[SwitchSSRN   ShearMod    PoissonRatio     R_Density  Crit_TooClose     TooCloseNormal_Multiplier MinNormalStress])
        write(io, "Ctr_X\tCtr_Y\tCtr_Z\tSt_L\tDip_L\tStAng\tDipAng\tLR/RN\ta\tb\tDc\tTheta_i\tV_i\tFric_i\tSig0\tSigGrad\tV_Const\tMaxLeng\n")
        writedlm(io, InputFile)
    end;

    println("Saved File Name: ",OutputFileName)
    end


BuildBulk()