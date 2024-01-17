
include("../Functions_BuildInputFile.jl")
using DelimitedFiles
using PyPlot
using PyCall
pygui(true)


function BuildBulk()


    OutputFileName="Input_BulkFaultGeometry.txt"

    SwitchSSRN = 1
    ShearMod= 20e9
    PoissonRatio = 0.2
    R_Density=2400.0
    Crit_TooClose= 1.05
    TooCloseNormal_Multiplier = 0.6
    MinNormalStress = 2e6
    

    Faultcount=3

    Fault1Center=[0.0,0.0,1500.0]
    Fault1StrikeAngle = 10.0
    Fault1DipAngle = 90.0
    Fault1Strikelength = 3000.0

    Fault2DipAngle = 90.0
    Fault2StrikeAngle = 160.0
    Fault2Strikelength = 1000.0
    Fault2Center=[Fault2Strikelength/2*cosd(20),-Fault2Strikelength/2*sind(20),1500]
    
    Theta_i = 1e3
    V_i = 1e-9

    MaxDiscritLength=200.0
    

    ############################# Write Bulk Input #################################
    ######++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++######
    ## Bulk File Order 
    ##  1.Ctr_X     2.Ctr_Y 3.Ctr_Z 4.St_L	    5.Dip_L	    6.StAng	    7.DipAng	8.LR/RN (-1: LL / +1 RL or -1:Reverse / +1 Normal)
    ##  9.a         10.b	11.Dc	12.Theta_i	13. V_i     14. Friction_i 15.NormalStress at surface [Pa]  
    ##  16. NoarmalStress Gradient [Pa] 17. V_Const     18. Minimum Segment Length



    ########################### Write Unstable Faults ##############################
    ######++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++######

    Fault = zeros(Faultcount,18)
    ### Fault 1.     
    Fault[1,1]=Fault1Center[1] + Fault1Strikelength/4 * cosd(Fault1StrikeAngle)
    Fault[1,2]=Fault1Center[2] + Fault1Strikelength/4 * sind(Fault1StrikeAngle)
    Fault[1,3]=Fault1Center[3] 
    Fault[2,1]=Fault1Center[1] - Fault1Strikelength/4 * cosd(Fault1StrikeAngle)
    Fault[2,2]=Fault1Center[2] - Fault1Strikelength/4 * sind(Fault1StrikeAngle)
    Fault[2,3]=Fault1Center[3] 
    Fault[1:2,4] .= Fault1Strikelength /2 
    Fault[1:2,5] .= Fault[1,3]*2;
    Fault[1:2,6] .=Fault1StrikeAngle;
    Fault[1:2,7] .=Fault1DipAngle;
    ### Fault 2     
    Fault[3,1:3]=Fault2Center
    Fault[3,4]=Fault2Strikelength
    Fault[3,5]=Fault[2,3]*2;
    Fault[3,6]=Fault2StrikeAngle;
    Fault[3,7]=Fault2DipAngle;

    ### Common Values 
    Fault[:,8] .= 1.0;
    Fault[:,9] .= 0.003;
    Fault[:,10] .= 0.006;
    Fault[:,11] .= 5e-4;
    Fault[:,12] .= Theta_i;
    Fault[:,13] .= V_i;
    Fault[:,14] .= 0.6
    Fault[:,15] .= 2e6
    Fault[:,16] .= 6000.0
    Fault[:,17] .= 0.0
    Fault[:,18] .= MaxDiscritLength

    LF_Distance = 1000.0
    LF_Depth = 7000.0
    LF_Length = 15000.0
    LoadingRate=1e-9
    Loading_Faults = zeros(2,18)

    Loading_Faults[1,1]=0.0;
    Loading_Faults[1,2]=LF_Distance;
    Loading_Faults[1,3]=LF_Depth/2;
    Loading_Faults[1,4]=LF_Length;
    Loading_Faults[1,5]=LF_Depth;

    Loading_Faults[2,1]=0.0;
    Loading_Faults[2,2]=-LF_Distance;
    Loading_Faults[2,3]=LF_Depth/2;
    Loading_Faults[2,4]=LF_Length;
    Loading_Faults[2,5]=LF_Depth

    ### Common Values 
    Loading_Faults[:,6].=0.0;
    Loading_Faults[:,7].=90.0;
    Loading_Faults[:,8].=-1;
    Loading_Faults[:,9].= 0.01;
    Loading_Faults[:,10].=0.005;
    Loading_Faults[:,11].=1e-3;
    Loading_Faults[:,12].=1e-3/LoadingRate;
    Loading_Faults[:,13].=LoadingRate;
    Loading_Faults[:,14].=0.6;
    Loading_Faults[:,15].=3e6;
    Loading_Faults[:,16].=0;
    Loading_Faults[:,17].=LoadingRate;
    Loading_Faults[:,18].=1e10; # This should be very large as we don't allow descritization for this



    FaultPlot_3D(Fault[:,1:3],Fault[:,4], Fault[:,5], Fault[:,6], Fault[:,7], Fault[:,8])
    xlabel("x")
    ylabel("y")
    # FaultPlot_3D(Stable_Faults[:,1:3],Stable_Faults[:,4], Stable_Faults[:,5], Stable_Faults[:,6], Stable_Faults[:,7], Stable_Faults[:,8])
    FaultPlot_3D(Loading_Faults[:,1:3],Loading_Faults[:,4], Loading_Faults[:,5], Loading_Faults[:,6], Loading_Faults[:,7], Loading_Faults[:,8])

    # InputFile=[Fault;Stable_Faults;Loading_Faults]
    InputFile=[Fault;Loading_Faults]

    FaultSegmentCount=zeros(1,2)
    for i in eachindex(InputFile[:,1])    
        FaultSegmentCount=[FaultSegmentCount;[ceil(InputFile[i,4]/InputFile[i,18]), ceil(InputFile[i,5]/InputFile[i,18])]']
    end
    FaultSegmentCount=FaultSegmentCount[2:end,:]
    TotalSegmentCount=sum(Int64,FaultSegmentCount[:,1].*FaultSegmentCount[:,2])
    println("Discretized Element Count will be ", TotalSegmentCount)
    

    # println("Total ",size(InputFile,1)," Faults Generated (", size(Fault,1), " Rupture zone faults, ",size(Stable_Faults,1), " Stable fault, ", size(Loading_Faults,1)," Loading Faults)" )

    println("Total ",size(InputFile,1)," Faults Generated (", size(Fault,1), " Rupture zone faults)" )

    open(OutputFileName, "w") do io
        write(io,"SwitchSS/RN\tShearMod\tPoissonRatio\tR_Density\tCrit_TooClose\tTooCloseNormal_Multiplier\tMinimum_NS\n")
        writedlm(io,[SwitchSSRN   ShearMod    PoissonRatio     R_Density  Crit_TooClose     TooCloseNormal_Multiplier MinNormalStress])
        write(io, "Ctr_X\tCtr_Y\tCtr_Z\tSt_L\tDip_L\tStAng\tDipAng\tLR/RN\ta\tb\tDc\tTheta_i\tV_i\tFric_i\tSig0\tSigGrad\tV_Const\tMaxLeng\n")
        writedlm(io, InputFile)
    end;

    println("Saved File Name: ",OutputFileName)
end


BuildBulk()