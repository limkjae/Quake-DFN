
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

    Faultlength = 5e3
    CurveThickness = 1000.0 
    FaultDepthCenter=2000.0

    NumberOfDivision=15
    BulkFaultcount=NumberOfDivision
    
    # YPointForCos = collect(0.0:pi*2/NumberOfDivision:pi*2)
    YPointForCos = collect(pi:pi/NumberOfDivision:pi*2)
    XPoint =  (cos.(YPointForCos) .+1) *CurveThickness /2
    YPoint = YPointForCos / 2 / pi * Faultlength .- Faultlength / 2
    
    
    FaultCenterY = zeros(NumberOfDivision+1)
    FaultCenterX = zeros(NumberOfDivision+1)
    FaultStrikeAngle = zeros(NumberOfDivision+1)
    FaultLength = zeros(NumberOfDivision+1)
    for i=1:NumberOfDivision
        FaultCenterY[i] = (YPoint[i]+YPoint[i+1])/2
        FaultCenterX[i] = (XPoint[i]+XPoint[i+1])/2
        FaultLength[i] = sqrt((YPoint[i+1]-YPoint[i])^2 + (XPoint[i+1]-XPoint[i])^2)
        FaultStrikeAngle[i] = atand((YPoint[i+1]-YPoint[i]) / (XPoint[i+1]-XPoint[i]))
        if FaultStrikeAngle[i] <0
            FaultStrikeAngle[i] = FaultStrikeAngle[i] + 180
        end
    end
    
    FaultCenterY[NumberOfDivision+1] = -Faultlength/4
    FaultCenterX[NumberOfDivision+1] = 0.0
    FaultLength[NumberOfDivision+1] = Faultlength/2
    FaultStrikeAngle[NumberOfDivision+1] = 90.0 



    # for FaultIdx in eachindex(X1)
    #     PyPlot.plot([X1[FaultIdx],X2[FaultIdx]],[Y1[FaultIdx],Y2[FaultIdx]])
    # end


    F_LR = 1.0
    F_a = 0.003 
    F_b =0.006 
    F_Dc = 3e-4 
    F_Theta = 1e3
    F_Vini = 1e-15
    NormalStressSurface=2e6 #10e6
    NormalStressGrad=7000.0 #Pascal per meter3
    MaximumSegLength=maximum(FaultLength[1:NumberOfDivision])

    ############################# Write Bulk Input #################################
    ######++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++######
    ### Bulk File Order 
    ###  1.Ctr_X     2.Ctr_Y 3.Ctr_Z 4.St_L	    5.Dip_L	    6.StAng	    7.DipAng	8.LR/RN (-1: LL / +1 RL or -1:Reverse / +1 Normal)
    ###  9.a         10.b	11.Dc	12.Theta_i	13. V_i     14. Friction_i 15.NormalStress at surface [Pa]  
    ###  16. NoarmalStress Gradient [Pa] 17. V_Const     18. Minimum Segment Length


    FaultInfo=zeros(BulkFaultcount+1,18)
    FaultInfo[:,1] = FaultCenterX
    FaultInfo[:,2] = FaultCenterY
    FaultInfo[:,3] .= FaultDepthCenter;
    FaultInfo[:,4] = FaultLength
    FaultInfo[:,5] .= FaultDepthCenter*2
    FaultInfo[:,6] = FaultStrikeAngle;
    FaultInfo[:,7] .= 90.0;
    FaultInfo[:,8] .= F_LR;
    FaultInfo[:,9] .= F_a;
    FaultInfo[:,10] .= F_b;
    FaultInfo[:,11] .= F_Dc;
    FaultInfo[:,12] .= F_Theta;
    FaultInfo[:,13] .= F_Vini;
    FaultInfo[:,14] .= 0.6;
    FaultInfo[:,15] .= NormalStressSurface;
    FaultInfo[:,16] .= NormalStressGrad;
    FaultInfo[:,17] .= 0.0;
    FaultInfo[:,18] .= MaximumSegLength;

    FaultPlot_3D(FaultInfo[:,1:3],FaultInfo[:,4], FaultInfo[:,5], FaultInfo[:,6], FaultInfo[:,7], FaultInfo[:,8])



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
