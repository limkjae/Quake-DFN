

include("../Functions_BuildInputFile.jl")
using DelimitedFiles
using PyPlot
using PyCall
pygui(true)

function BuildBulk()


    OutputFileName="Input_BulkFaultGeometry.txt"
    
    SwitchSSRN = 1
    ShearMod= 32.038e9
    PoissonRatio = 0.25
    R_Density=2670.0
    Crit_TooClose= 1.05
    TooCloseNormal_Multiplier = 0.6
    MinNormalStress = 2e6

    VW_Thickness=12e3
    TZ_Length = 2e3
    RZ_HalfL=  32e3 # half length of Rupture zone
    RZ_Depth = TZ_Length+TZ_Length+TZ_Length+ VW_Thickness# depth of Rupture zone
    FZ_HalfL = 50e3 # half length of Fault zone (seismic + aseismic)
    FZ_Depth = 40e3 # Depth of Fault zone
    LF_HalfL = 1000e3 # Half length of Loading Fault
    LF_Depth = 1000e3 # depth of loading fault
    LF_Distance = 20e3 # Distance of loading faults from fault zone
    LoadingRate=1e-9
    NormalStressSurface=25e6 #10e6
    NormalStressGrad=0.0 #Pascal per meter


    Unstable_a=0.004
    Unstable_b=0.03
    Unstable_Dc=0.14
    Unstable_Theta=1e9
    UntableMaximumSegLength = 4000.0

    Stable_a=0.04
    Stable_b=0.03
    Stable_Dc=0.14 
    StableMaximumSegLength = 4000.0;


    ############################# Write Bulk Input #################################
    ######++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++######
    ## Bulk File Order 
    ##  1.Ctr_X     2.Ctr_Y 3.Ctr_Z 4.St_L	    5.Dip_L	    6.StAng	    7.DipAng	8.LR/RN (-1: LL / +1 RL or -1:Reverse / +1 Normal)
    ##  9.a         10.b	11.Dc	12.Theta_i	13. V_i     14. Friction_i 15.NormalStress at surface [Pa]  
    ##  16. NoarmalStress Gradient [Pa] 17. V_Const     18. Minimum Segment Length



    ########################### Write Unstable Faults ##############################
    ######++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++######

    # ### Unstable Fault 1. 
    Unstable_Faults=zeros(2,18)

    Unstable_Faults[1,1]=0;
    Unstable_Faults[1,2]=-24e3
    Unstable_Faults[1,3]=TZ_Length*2 + VW_Thickness/2;
    Unstable_Faults[1,4]=12e3;
    Unstable_Faults[1,5]=VW_Thickness

    Unstable_Faults[2,1]=0;
    Unstable_Faults[2,2]= 48e3/2 - 18e3
    Unstable_Faults[2,3]=TZ_Length*2 + VW_Thickness/2;
    Unstable_Faults[2,4]=48e3;
    Unstable_Faults[2,5]=VW_Thickness

    ### Common Values 
    Unstable_Faults[:,6].=90.0;
    Unstable_Faults[:,7].=90.0;
    Unstable_Faults[:,8].=1.0;
    Unstable_Faults[:,9].=Unstable_a;
    Unstable_Faults[:,10].=Unstable_b;
    Unstable_Faults[1,11]=0.13
    Unstable_Faults[2,11]=0.14
    Unstable_Faults[1,12]=0.13 / LoadingRate;
    Unstable_Faults[2,12]=0.14 / LoadingRate;
    Unstable_Faults[1,13]=3e-2;
    Unstable_Faults[2,13]=LoadingRate;
    Unstable_Faults[:,14].=0.6;
    Unstable_Faults[:,15].=NormalStressSurface;
    Unstable_Faults[:,16].=NormalStressGrad;
    Unstable_Faults[:,17].=0.0;
    Unstable_Faults[:,18].=UntableMaximumSegLength;

    
    ########################### Write Transition Faults 1 ############################
    ######++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++######
    ### Stable Fault 1. 
    Transition_Faults=zeros(8,18)

    Transition_Faults[1,1]=0;
    Transition_Faults[1,2]=RZ_HalfL - TZ_Length/4
    Transition_Faults[1,3]=TZ_Length + (RZ_Depth - TZ_Length)/2;
    Transition_Faults[1,4]=TZ_Length/2;
    Transition_Faults[1,5]=(RZ_Depth - TZ_Length)

    # ### Stable Fault 2. 

    Transition_Faults[2,1]=0;
    Transition_Faults[2,2]=-(RZ_HalfL - TZ_Length/4 )
    Transition_Faults[2,3]=TZ_Length + (RZ_Depth - TZ_Length)/2;
    Transition_Faults[2,4]=TZ_Length/2;
    Transition_Faults[2,5]=(RZ_Depth - TZ_Length)

    ### Stable Fault 3. 
    Transition_Faults[3,1]=0;
    Transition_Faults[3,2]=0;
    Transition_Faults[3,3]=TZ_Length+TZ_Length/4;
    Transition_Faults[3,4]=(RZ_HalfL-TZ_Length/2)*2;
    Transition_Faults[3,5]=TZ_Length/2;

    ### Stable Fault 4. 
    Transition_Faults[4,1]=0;
    Transition_Faults[4,2]=0;
    Transition_Faults[4,3]=RZ_Depth-TZ_Length/4;
    Transition_Faults[4,4]=(RZ_HalfL-TZ_Length/2)*2;
    Transition_Faults[4,5]=TZ_Length/2;

    ### Stable Fault 5. 
    Transition_Faults[5,1]=0;
    Transition_Faults[5,2]=RZ_HalfL - TZ_Length*3/4 
    Transition_Faults[5,3]=TZ_Length + (RZ_Depth - TZ_Length)/2;
    Transition_Faults[5,4]=TZ_Length/2;
    Transition_Faults[5,5]=(RZ_Depth - TZ_Length*2)

    # ### Stable Fault 6. 

    Transition_Faults[6,1]=0;
    Transition_Faults[6,2]=-(RZ_HalfL - TZ_Length*3/4 )
    Transition_Faults[6,3]=TZ_Length + (RZ_Depth - TZ_Length)/2;
    Transition_Faults[6,4]=TZ_Length/2;
    Transition_Faults[6,5]=(RZ_Depth - TZ_Length*2)

    ### Stable Fault 7. 
    Transition_Faults[7,1]=0;
    Transition_Faults[7,2]=0;
    Transition_Faults[7,3]=TZ_Length+TZ_Length*3/4;
    Transition_Faults[7,4]=(RZ_HalfL-TZ_Length)*2;
    Transition_Faults[7,5]=TZ_Length/2;

    ### Stable Fault 8. 
    Transition_Faults[8,1]=0;
    Transition_Faults[8,2]=0;
    Transition_Faults[8,3]=RZ_Depth-TZ_Length*3/4;
    Transition_Faults[8,4]=(RZ_HalfL-TZ_Length)*2;
    Transition_Faults[8,5]=TZ_Length/2;


    ### Common Values 
    Transition_Faults[:,6].=90.0;
    Transition_Faults[:,7].=90.0;
    Transition_Faults[:,8].=1.0;
    Transition_Faults[1:4,9].=Unstable_a + (Stable_a - Unstable_a)*3/4;
    Transition_Faults[5:8,9].=Unstable_a + (Stable_a - Unstable_a)*1/4;
    Transition_Faults[:,10].=Stable_b;
    Transition_Faults[:,11].=Stable_Dc;
    Transition_Faults[:,12].=Stable_Dc/LoadingRate;
    Transition_Faults[:,13].=LoadingRate;
    Transition_Faults[:,14].=0.6;
    Transition_Faults[:,15].=NormalStressSurface;
    Transition_Faults[:,16].=NormalStressGrad;
    Transition_Faults[:,17].=0.0;
    Transition_Faults[:,18].=UntableMaximumSegLength;



    ###################### Write Stable Surrounding Faults #########################
    ######++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++######
    ### Stable Fault 1. 
    Stable_Faults=zeros(4,18)
    Stable_Faults[1,1]=0;
    Stable_Faults[1,2]=RZ_HalfL+(FZ_HalfL-RZ_HalfL)/2;
    Stable_Faults[1,3]=RZ_Depth/2;
    Stable_Faults[1,4]=(FZ_HalfL-RZ_HalfL);
    Stable_Faults[1,5]=RZ_Depth;
    ### Stable Fault 2. 
    Stable_Faults[2,1]=0;
    Stable_Faults[2,2]=-(RZ_HalfL+(FZ_HalfL-RZ_HalfL)/2)
    Stable_Faults[2,3]=RZ_Depth/2;
    Stable_Faults[2,4]=(FZ_HalfL-RZ_HalfL);
    Stable_Faults[2,5]=RZ_Depth;
    ### Stable Fault 3. 
    Stable_Faults[3,1]=0;
    Stable_Faults[3,2]=0;
    Stable_Faults[3,3]=RZ_Depth+(FZ_Depth-RZ_Depth)/2;
    Stable_Faults[3,4]=FZ_HalfL*2;
    Stable_Faults[3,5]=FZ_Depth-RZ_Depth;

    ### Stable Fault 4. 
    Stable_Faults[4,1]=0;
    Stable_Faults[4,2]=0;
    Stable_Faults[4,3]=TZ_Length/2;
    Stable_Faults[4,4]=RZ_HalfL*2;
    Stable_Faults[4,5]=TZ_Length;


    ### Common Values 
    Stable_Faults[:,6].=90.0;
    Stable_Faults[:,7].=90.0;
    Stable_Faults[:,8].=1.0;
    Stable_Faults[:,9].=Stable_a;
    Stable_Faults[:,10].=Stable_b;
    Stable_Faults[:,11].=Stable_Dc;
    Stable_Faults[:,12].=Stable_Dc/LoadingRate;
    Stable_Faults[:,13].=LoadingRate;
    Stable_Faults[:,14].=0.6;
    Stable_Faults[:,15].=NormalStressSurface;
    Stable_Faults[:,16].=NormalStressGrad;
    Stable_Faults[:,17].=0.0;
    Stable_Faults[1:3,18].=StableMaximumSegLength;
    Stable_Faults[4,18]=StableMaximumSegLength;


    ##  1.Ctr_X     2.Ctr_Y 3.Ctr_Z 4.St_L	    5.Dip_L	    6.StAng	    7.DipAng	8.LR
    ##  9.a         10.b	11.Dc	12.Theta_i	13. V_i     14. Friction_i 15.NormalStress at surface [Pa]  
    ##  16. NoarmalStress Gradient [Pa] 17. V_Const     18. Minimum Segment Length


    Loading_Faults=zeros(5,18)

    Loading_Faults[1,1]=0;
    Loading_Faults[1,2]=FZ_HalfL+(LF_HalfL-FZ_HalfL)/2;
    Loading_Faults[1,3]=FZ_Depth/2;
    Loading_Faults[1,4]=(LF_HalfL-FZ_HalfL);
    Loading_Faults[1,5]=FZ_Depth;

    Loading_Faults[2,1]=0;
    Loading_Faults[2,2]=-(FZ_HalfL+(LF_HalfL-FZ_HalfL)/2)
    Loading_Faults[2,3]=FZ_Depth/2;
    Loading_Faults[2,4]=(LF_HalfL-FZ_HalfL);
    Loading_Faults[2,5]=FZ_Depth;

    Loading_Faults[3,1]=0;
    Loading_Faults[3,2]=0;
    Loading_Faults[3,3]=FZ_Depth+(LF_Depth-FZ_Depth)/2;
    Loading_Faults[3,4]=LF_HalfL*2;
    Loading_Faults[3,5]=LF_Depth-FZ_Depth;

    Loading_Faults[4,1]=LF_Distance;
    Loading_Faults[4,2]=0;
    Loading_Faults[4,3]=LF_Depth/2;
    Loading_Faults[4,4]=LF_HalfL*2;
    Loading_Faults[4,5]=LF_Depth;

    Loading_Faults[5,1]=-LF_Distance;
    Loading_Faults[5,2]=0;
    Loading_Faults[5,3]=LF_Depth/2;
    Loading_Faults[5,4]=LF_HalfL*2;
    Loading_Faults[5,5]=LF_Depth;

    ### Common Values 
    Loading_Faults[:,6].=90.0;
    Loading_Faults[:,7].=90.0;
    Loading_Faults[1:3,8].=1.0;
    Loading_Faults[4:5,8].=-1.0;
    Loading_Faults[:,9].=Stable_a;
    Loading_Faults[:,10].=Stable_b;
    Loading_Faults[:,11].=Stable_Dc;
    Loading_Faults[:,12].=Stable_Dc/LoadingRate;
    Loading_Faults[:,13].=LoadingRate;
    Loading_Faults[:,14].=0.6;
    Loading_Faults[:,15].=NormalStressSurface;
    Loading_Faults[:,16].=NormalStressGrad;
    Loading_Faults[:,17].=LoadingRate;
    Loading_Faults[:,18].=1e10; # This should be very large as we don't allow descritization for this


    FaultPlot_3D(Unstable_Faults[:,1:3],Unstable_Faults[:,4], Unstable_Faults[:,5], Unstable_Faults[:,6], Unstable_Faults[:,7], Unstable_Faults[:,8])
    FaultPlot_3D(Transition_Faults[:,1:3],Transition_Faults[:,4], Transition_Faults[:,5], Transition_Faults[:,6], Transition_Faults[:,7], Transition_Faults[:,8])
    FaultPlot_3D(Stable_Faults[:,1:3],Stable_Faults[:,4], Stable_Faults[:,5], Stable_Faults[:,6], Stable_Faults[:,7], Stable_Faults[:,8])
    FaultPlot_3D(Loading_Faults[:,1:3],Loading_Faults[:,4], Loading_Faults[:,5], Loading_Faults[:,6], Loading_Faults[:,7], Loading_Faults[:,8])

    # InputFile=[Unstable_Faults;Stable_Faults;Loading_Faults]
    InputFile=[Unstable_Faults; Transition_Faults; Stable_Faults; Loading_Faults]

    # FaultPlot_3D(InputFile[:,1:3],InputFile[:,4], InputFile[:,5], InputFile[:,6], InputFile[:,7], InputFile[:,8], InputFile[:,9]-InputFile[:,10])

    # FaultPlot_3D(Unstable_Faults[:,1:3],Unstable_Faults[:,4], Unstable_Faults[:,5], Unstable_Faults[:,6], Unstable_Faults[:,7], Unstable_Faults[:,8])


    FaultSegmentCount=zeros(1,2)
    for i in eachindex(InputFile[:,1])    
        FaultSegmentCount=[FaultSegmentCount;[ceil(InputFile[i,4]/InputFile[i,18]), ceil(InputFile[i,5]/InputFile[i,18])]']
    end
    FaultSegmentCount=FaultSegmentCount[2:end,:]
    TotalSegmentCount=sum(Int64,FaultSegmentCount[:,1].*FaultSegmentCount[:,2])
    println("Discretized Element Count will be ", TotalSegmentCount)
    

    println("Total ",size(InputFile,1)," Faults Generated (", size(Stable_Faults,1)+size(Transition_Faults,1)+size(Unstable_Faults,1), " Regular faults, ", size(Loading_Faults,1)," Loading Faults)" )

    open(OutputFileName, "w") do io
        write(io,"SwitchSS/RN\tShearMod\tPoissonRatio\tR_Density\tCrit_TooClose\tTooCloseNormal_Multiplier\tMinimum_NS\n")
        writedlm(io,[SwitchSSRN   ShearMod    PoissonRatio     R_Density  Crit_TooClose     TooCloseNormal_Multiplier MinNormalStress])
        write(io, "Ctr_X\tCtr_Y\tCtr_Z\tSt_L\tDip_L\tStAng\tDipAng\tLR/RN\ta\tb\tDc\tTheta_i\tV_i\tFric_i\tSig0\tSigGrad\tV_Const\tMaxLeng\n")
        writedlm(io, InputFile)
    end;

    println("Saved File Name: ",OutputFileName)
end


BuildBulk()

#=
TotalStep=300000; # Total simulation step
RecordStep=100; # Simulation sampling rate 
SwitchV=3e-5; # velocity criteria where slow and fast solution method switches
TimeStepping =
    [300000.0 500e-1 100e-1 0.5e-1;
    1e-8   3e-5  1e-3  1e-2]


ShearModulus=3464.0^2*2670.0
PoissonRatio=0.25
RockDensity=2670.0

=#