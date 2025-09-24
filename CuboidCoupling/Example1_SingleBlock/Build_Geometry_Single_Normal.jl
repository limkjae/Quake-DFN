include("../../Functions_BuildInputFile.jl")
using DelimitedFiles
using JLD2
using PyPlot
using PyCall
pygui(false)

 
function BuildBulk()
    OutputFileName="Input_BulkFaultGeometry.txt"

    # Fault parameters
    Phi = 70 # deg (Conventional Dip Angle)
    DipAngle = 180.0 - Phi # deg (QuakeDFN Dip Angle)
    StrikeAngle = 90.0 # deg
    FaultCenterX = 0
    FaultCenterY = 0
    FaultCenterZ = 3000
    LengthDip = 800
    LengthStrike = 2000

    println("Dip Angle in QuakeDFN is ", DipAngle)
    println("Fault Center is: ",FaultCenterX, "  ",FaultCenterY, "  ", FaultCenterZ )
    
    MyPossionRatio = 0.15
    MyYoungModulus = 15e9
    MyShearModulus = MyYoungModulus/(2*MyPossionRatio+2)

    SwitchSSRN = 2 # 3D Normal fault
    ShearMod = MyShearModulus
    PoissonRatio = MyPossionRatio
    R_Density = 2400.0

    Crit_TooClose= 1.05
    TooCloseNormal_Multiplier = 0.6
    
    MinNormalStress = 2e6
    NormalStressSurface = 0.0   #18e6
    NormalStressGrad = 0.0       # 7000.0 #Pascal per meter

    UnstableFaultcount = 1

    Unstable_a=0.003
    Unstable_b=0.006
    Unstable_Dc= 1e-4
    Unstable_Theta=1e3
    UnstableMaximumSegLength = 20.0 # 200.0 | 50.0 | 20.0


    ############################# Write Bulk Input #################################
    ######++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++######
    ## Bulk File Order 
    ##  1.Ctr_X     2.Ctr_Y 3.Ctr_Z 4.St_L	    5.Dip_L	    6.StAng	    7.DipAng	8.LR/RN (-1: LL / +1 RL or -1:Reverse / +1 Normal)
    ##  9.a         10.b	11.Dc	12.Theta_i	13. V_i     14. Friction_i 15.NormalStress at surface [Pa]  
    ##  16. NormalStress Gradient [Pa] 17. V_Const     18. Minimum Segment Length



    ########################### Write Unstable Faults ##############################
    ######++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++######

    ### Stable Fault 1. 

    Unstable_Faults=zeros(UnstableFaultcount,18)
    Unstable_Faults[:,1].=FaultCenterX;
    Unstable_Faults[:,2].=FaultCenterY;
    Unstable_Faults[:,3].=FaultCenterZ;
    Unstable_Faults[:,4].=LengthStrike; 
    Unstable_Faults[:,5].=LengthDip;
    Unstable_Faults[:,6].=StrikeAngle;
    Unstable_Faults[:,7].=DipAngle;

    ### Common Values 
    Unstable_Faults[:,8].=1.0;
    Unstable_Faults[:,9].=Unstable_a;
    Unstable_Faults[:,10].=Unstable_b;
    Unstable_Faults[:,11].=Unstable_Dc;
    Unstable_Faults[:,12].=Unstable_Theta;
    Unstable_Faults[:,13].=1e-10;
    Unstable_Faults[:,14].=0.6; # 0.3
    Unstable_Faults[:,15].=NormalStressSurface;
    Unstable_Faults[:,16].=NormalStressGrad;
    Unstable_Faults[:,17].=0.0;
    Unstable_Faults[:,18].=UnstableMaximumSegLength;

    FaultPlot_3D(Unstable_Faults[:,1:3],Unstable_Faults[:,4], Unstable_Faults[:,5], Unstable_Faults[:,6], Unstable_Faults[:,7], Unstable_Faults[:,8])
    xlabel("x (m)")
    ylabel("y (m)")

    InputFile=[Unstable_Faults;]
    FaultSegmentCount=zeros(1,2)
    for i in eachindex(InputFile[:,1])    
        FaultSegmentCount=[FaultSegmentCount;[ceil(InputFile[i,4]/InputFile[i,18]), ceil(InputFile[i,5]/InputFile[i,18])]']

    end
    FaultSegmentCount=FaultSegmentCount[2:end,:]

    TotalSegmentCount=sum(Int64,FaultSegmentCount[:,1].*FaultSegmentCount[:,2])
    println("Discretized Element Count will be ", TotalSegmentCount)


    println("Total ",size(InputFile,1)," Faults Generated")

    open(OutputFileName, "w") do io
        write(io,"SwitchSS/RN\tShearMod\tPoissonRatio\tR_Density\tCrit_TooClose\tTooCloseNormal_Multiplier\tMinimum_NS\n")
        writedlm(io,[SwitchSSRN   ShearMod    PoissonRatio     R_Density  Crit_TooClose     TooCloseNormal_Multiplier MinNormalStress])
        write(io, "Ctr_X\tCtr_Y\tCtr_Z\tSt_L\tDip_L\tStAng\tDipAng\tLR/RN\ta\tb\tDc\tTheta_i\tV_i\tFric_i\tSig0\tSigGrad\tV_Const\tMaxLeng\n")
        writedlm(io, InputFile)
    end;

    println("Saved File Name: ",OutputFileName)
end

BuildBulk()
