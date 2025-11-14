
include("../../scripts/Functions_BuildInputFile.jl")
using DelimitedFiles
using PyPlot
using PyCall
pygui(true)


#= ######################################################################
This simulation is for benchmarking against the step over simulation 
Conducted by Kroll et al., (2023). More details can be found here:

Kayla A Kroll, James H Dieterich, Keith B Richards-Dinger, David D Oglesby, 
3-D Simulations of earthquakes rupture jumps: 1. Homogeneous pre-stress conditions, 
Geophysical Journal International, Volume 234, Issue 1, July 2023, Pages 395–403,
https://doi.org/10.1093/gji/ggad048
####################################################################### =#


#= ####################################################################################################
##### Following  Lines should be included in the QuickParameterAdjust.jl to initiate the rupture ######

NucleationPoint=[-10e3,0, 10e3]
Friction0 = Fault_Friction_i[1] - Fault_a[1] * log.(Fault_V_i[1]/1e-9) - Fault_b[1] * log(Fault_Theta_i[1] * 1e-9/Fault_Dc[1]); # Initial friction
# println(Friction0)
for i=1:FaultCount        if norm(NucleationPoint - FaultCenter[i,:]) < 3e3
        Fault_V_i[i] = 1e-6
        Fault_Theta_i[i] = exp(((Fault_Friction_i[i] - Friction0) - Fault_a[i] * log(Fault_V_i[i]/1e-9))/Fault_b[i])/1e-9*Fault_Dc[i]
     end
end
###############-----------------------------------------------------------------------################
################################################################################################### =#

function BuildBulk()


    OutputFileName="Input_BulkFaultGeometry.txt"

    SwitchSSRN = 1
    R_Density=2670.0
    ShearMod= 32e9
    PoissonRatio = 0.25
    Crit_TooClose= 1.05
    TooCloseNormal_Multiplier = 0.6
    MinNormalStress = 2e6
    


    Fault1StrikeAngle = 0.0
    Fault1DipAngle = 90.0

    StepOver = -0.6e3
    
    Faultcount=2
    MaxDiscritLength=1000.0
    
    RSF_a = 0.01
    RSF_b = 0.012
    RSF_Dc = 1e-5
    Theta_i = 2.6e10
    InitialShearStress = 29.38e6    
    MaxSigSurface = 60e6
    MaxSigGrad = 0.0 #Pa/m
    V_i = 2.17e-13
    Mu_i = InitialShearStress/MaxSigSurface
    

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
    Fault[1,1]=-10000.0
    Fault[1,2]=-0.0;
    Fault[1,3]=20000/2;
    Fault[1,4]=30e3;
    Fault[1,5]=Fault[1,3]*2;
    Fault[1,6]=Fault1StrikeAngle;
    Fault[1,7]=Fault1DipAngle;

    ### Fault 2     
    Fault[2,1]=10000.0
    Fault[2,2]=StepOver;
    Fault[2,3]=20000/2;
    Fault[2,4]=30e3
    Fault[2,5]=Fault[1,3]*2;
    Fault[2,6]=Fault1StrikeAngle;
    Fault[2,7]=Fault1DipAngle;
    
    ### Common Values 
    Fault[:,9] .= RSF_a;
    Fault[:,10] .= RSF_b;
    Fault[:,11] .= RSF_Dc;
    Fault[:,12] .= Theta_i;
    Fault[:,13] .= V_i;

    ### Friction and NormalSTress
    
    ############################# Write Bulk Input #################################
    ######++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++######
    ## Bulk File Order 
    ##  1.Ctr_X     2.Ctr_Y 3.Ctr_Z 4.St_L	    5.Dip_L	    6.StAng	    7.DipAng	8.LR/RN (-1: LL / +1 RL or -1:Reverse / +1 Normal)
    ##  9.a         10.b	11.Dc	12.Theta_i	13. V_i     14. Friction_i 15.NormalStress at surface [Pa]  
    ##  16. NoarmalStress Gradient [Pa] 17. V_Const     18. Minimum Segment Length

    Fault[:,8].= -1;
    Fault[:,14] .= Mu_i;
    Fault[:,15] .= MaxSigSurface #.* NormalStressParameter;
    Fault[:,16] .= MaxSigGrad;
    Fault[:,17] .= 0.0;
    Fault[:,18] .= MaxDiscritLength;
    


    
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