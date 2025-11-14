
include("../../scripts/Functions_BuildInputFile.jl")
using DelimitedFiles
using PyPlot
using PyCall
pygui(true)


function BuildBulk()


    OutputFileName="Input_BulkFaultGeometry.txt"

    SwitchSSRN = 2
    ShearMod= 20e9
    PoissonRatio = 0.2
    R_Density=2400.0
    Crit_TooClose= 1.05
    TooCloseNormal_Multiplier = 0.6
    MinNormalStress = 2e6
    

    NormalStressSurface=0 #10e6

    Fault1StrikeAngle = 90.0
    Fault2StrikeAngle = 90.0
    Fault1DipAngle = 60.0
    Fault2DipAngle = 110.0
    MinFriction = 0.2
    
    Faultcount=3
    MaxDiscritLength=200.0
    
    MaxSigSurface = 2e6
    MaxSigGrad = 6000.0 #Pa/m

    MaxStressOrientation = 45. # between 0-180 degree
    StressRatioMaxOverMin = 0.5
    ############################# Write Bulk Input #################################
    ######++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++######
    ## Bulk File Order 
    ##  1.Ctr_X     2.Ctr_Y 3.Ctr_Z 4.St_L	    5.Dip_L	    6.StAng	    7.DipAng	8.LR/RN (-1: LL / +1 RL or -1:Reverse / +1 Normal)
    ##  9.a         10.b	11.Dc	12.Theta_i	13. V_i     14. Friction_i 15.NormalStress at surface [Pa]  
    ##  16. NoarmalStress Gradient [Pa] 17. V_Const     18. Minimum Segment Length



    ########################### Write Unstable Faults ##############################
    ######++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++######

    Fault=zeros(Faultcount,18)
    ### Fault 1.     
    Fault[1,1] = - 2000 / tand(Fault1DipAngle) / 2
    Fault[1,2]=-0.0;
    Fault[1,3]=1000.;
    Fault[1,4]= 2000;
    Fault[1,5]= 2000 / sind(Fault1DipAngle) ;
    Fault[1,6]=Fault1StrikeAngle;
    Fault[1,7]=Fault1DipAngle;

    ### Fault 2     
    Fault[2,1]= 2000 / tand(180 - Fault2DipAngle) / 2
    Fault[2,2]=0.0;
    Fault[2,3]=1000.;
    Fault[2,4]=2000.;
    Fault[2,5]=2000 / sind(180 - Fault2DipAngle) 
    Fault[2,6]=Fault2StrikeAngle;
    Fault[2,7]=Fault2DipAngle;

    ### Fault 3     
    Fault[3,1]= -2000 / tand(180 - Fault2DipAngle) / 2
    Fault[3,2]=0.0;
    Fault[3,3]=3000.;
    Fault[3,4]=2000.;
    Fault[3,5]=2000 / sind(180 - Fault2DipAngle) 
    Fault[3,6]=Fault2StrikeAngle;
    Fault[3,7]=Fault2DipAngle;

    Fault[:,8].= -1
    ### Common Values 
    Fault[:,9] .= 0.003;
    Fault[:,10] .= 0.006;
    Fault[:,11] .= 2e-4;
    Fault[:,12] .= 1e10;
    Fault[:,13] .= 1e-15;
    
    Fault[:,14] .= 0.33;
    Fault[:,15] .= MaxSigSurface #.* NormalStressParameter;
    Fault[:,16] .= 4500.0;
    Fault[:,17].=0.0;
    Fault[:,18].=MaxDiscritLength;

    
    FaultPlot_3D(Fault[:,1:3],Fault[:,4], Fault[:,5], Fault[:,6], Fault[:,7], Fault[:,8])
    xlabel("x")
    ylabel("y")

    InputFile=[Fault;]

    FaultSegmentCount=zeros(1,2)
    for i in eachindex(InputFile[:,1])    
        FaultSegmentCount=[FaultSegmentCount;[ceil(InputFile[i,4]/InputFile[i,18]), ceil(InputFile[i,5]/InputFile[i,18])]']
    end
    FaultSegmentCount=FaultSegmentCount[2:end,:]
    TotalSegmentCount=sum(Int64,FaultSegmentCount[:,1].*FaultSegmentCount[:,2])
    println("TotalSegment Count will be ", TotalSegmentCount)
    

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