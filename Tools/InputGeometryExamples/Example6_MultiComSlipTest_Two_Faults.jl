
include("../../scripts/Functions_BuildInputFile.jl")
using DelimitedFiles
using PyPlot
using PyCall
pygui(true)


function BuildBulk()


    OutputFileName="Input_BulkFaultGeometry.txt"

    SwitchSSRN = 0
    ShearMod= 20e9
    PoissonRatio = 0.2
    R_Density=2400.0
    Crit_TooClose= 1.05
    TooCloseNormal_Multiplier = 0.6
    MinNormalStress = 2e6
    

    NormalStressSurface=0 #10e6

    Fault1StrikeAngle = 30.0
    Fault2StrikeAngle = 110.0
    Fault1DipAngle = 70.0
    Fault2DipAngle = 45.0
    Fault1RakeAngle = 90.0
    Fault2RakeAngle = 30.0

    MinFriction = 0.1
    
    Faultcount=2
    MaxDiscritLength=200.0
    
    MaxSigSurface = 2e6
    MaxSigGrad = 6000.0 #Pa/m

    MaxStressOrientation = 45. # between 0-180 degree
    StressRatioMaxOverMin = 0.5
    ############################# Write Bulk Input #################################
    ######++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++######
    ## Bulk File Order 
    ##  1.Ctr_X     2.Ctr_Y 3.Ctr_Z 4.St_L	    5.Dip_L	    6.StAng	    7.DipAng	8.Rake
    ##  9.a         10.b	11.Dc	12.Theta_i	13. V_i     14. Friction_i 15.NormalStress at surface [Pa]  
    ##  16. NoarmalStress Gradient [Pa] 17. V_Const     18. Minimum Segment Length



    ########################### Write Unstable Faults ##############################
    ######++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++######

    Fault=zeros(Faultcount,18)
    ### Fault 1.     
    Fault[1,1]=-1000.0
    Fault[1,2]=-0.0;
    Fault[1,3]=1500.;
    Fault[1,4]=2000.;
    Fault[1,5]=Fault[1,3]*2;
    Fault[1,6]=Fault1StrikeAngle;
    Fault[1,7]=Fault1DipAngle;
    Fault[1,8]=Fault1RakeAngle;

    ### Fault 2     
    Fault[2,1]=100.0
    Fault[2,2]=000.0;
    Fault[2,3]=1500.;
    Fault[2,4]=4000. / 2;
    Fault[2,5]=Fault[2,3]*2;
    Fault[2,6]=Fault2StrikeAngle;
    Fault[2,7]=Fault2DipAngle;
    Fault[2,8]=Fault2RakeAngle;

    # ### Fault 3     
    # Fault[3,1]=0.0
    # Fault[3,2]=-1000.0;
    # Fault[3,3]=1500.;
    # Fault[3,4]=4000. / 2;
    # Fault[3,5]=Fault[3,3]*2;
    # Fault[3,6]=Fault2StrikeAngle;
    # Fault[3,7]=Fault2DipAngle;

    ### Common Values 
    Fault[:,9] .= 0.003;
    Fault[:,10] .= 0.006;
    Fault[:,11] .= 2e-4;
    Fault[:,12] .= 1e10;
    Fault[:,13] .= 1e-15;

    ### Friction and NormalSTress
    # NormalStressParameter = (1+StressRatioMaxOverMin)/2 .+ (1-StressRatioMaxOverMin)/2 .* cosd.(2 .* (Fault[:,6] .- 90.0 .- MaxStressOrientation))
    # ShearStressParameter = -(1-StressRatioMaxOverMin)/2 .* sind.(2 .* (Fault[:,6] .- 90.0 .- MaxStressOrientation))
    # Fault[:,14] = abs.(ShearStressParameter ./ NormalStressParameter);
    Fault[:,14] .= 0.6;
    
    Fault[:,15] .= MaxSigSurface #.* NormalStressParameter;
    Fault[:,16] .= 6000.0;
    Fault[:,17].=0.0;
    Fault[:,18].=MaxDiscritLength;
   
    
    for i=1:length(Fault[:,14])
        if Fault[i,14] < MinFriction
            Fault[i,14] = MinFriction
        end
    end
    # println(ShearStressParameter,"  ", NormalStressParameter, "  ", abs.(ShearStressParameter ./ NormalStressParameter))
    # Fault1NS = (1+StressRatioMaxOverMin)/2 + (1-StressRatioMaxOverMin)/2*cosd(2*(Fault1StrikeAngle-90-MaxStressOrientation))
    # Fault1SS = (1-StressRatioMaxOverMin)/2*sind(2*(Fault1StrikeAngle-90-MaxStressOrientation))
    
    


    figure(3)
    clf()
    FaultPlot_3D(Fault[:,1:3],Fault[:,4], Fault[:,5], Fault[:,6], Fault[:,7], Fault[:,8])
    xlabel("x")
    ylabel("y")
    # FaultPlot_3D(Stable_Faults[:,1:3],Stable_Faults[:,4], Stable_Faults[:,5], Stable_Faults[:,6], Stable_Faults[:,7], Stable_Faults[:,8])
    # FaultPlot_3D(Loading_Faults[:,1:3],Loading_Faults[:,4], Loading_Faults[:,5], Loading_Faults[:,6], Loading_Faults[:,7], Loading_Faults[:,8])

    # InputFile=[Fault;Stable_Faults;Loading_Faults]
    InputFile=[Fault;]

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
        write(io, "Ctr_X\tCtr_Y\tCtr_Z\tSt_L\tDip_L\tStAng\tDipAng\tRake\ta\tb\tDc\tTheta_i\tV_i\tFric_i\tSig0\tSigGrad\tV_Const\tMaxLeng\n")
        writedlm(io, InputFile)
    end;

    println("Saved File Name: ",OutputFileName)
end


BuildBulk()