
include("../Functions_BuildInputFile.jl")
using DelimitedFiles
using PyPlot
using PyCall
pygui(true)

#######################################################################
### This simulation is for benchmarking against the TPV18 of 
### The SCEC/USGS Spontaneous Rupture Code Verification Project
### See https://strike.scec.org/cvws/cgi-bin/cvws.cgi for more detail
#######################################################################


function BuildBulk()


    OutputFileName="Input_BulkFaultGeometry.txt"

    SwitchSSRN = 1
    R_Density=2700.0
    ShearMod= 3300^2 * R_Density
    BulkMod = 5716^2 * R_Density -4/3*ShearMod
    PoissonRatio = (3*BulkMod - 2*ShearMod)/2/(3*BulkMod + ShearMod)
    Crit_TooClose= 1.05
    TooCloseNormal_Multiplier = 0.6
    MinNormalStress = 2e6
    
    
    Fault1StrikeAngle = 0.0
    Fault2StrikeAngle = 150.0
    Fault1DipAngle = 90.0
    Fault2DipAngle = 90.0
    
    Faultcount=7
    MaxDiscritLength=1000.0
    
    RSF_a = 0.006
    RSF_b = 0.012
    RSF_Dc = 0.001

    NormalStressI =2e6
    NormalStressGradV = 2700*9.8
    PressureGrad = 1000*9.8
    RotationAngle = -30
    StressGradTensor = 
    [0.50911055*(NormalStressGradV-PressureGrad)  -0.15487679* (NormalStressGradV-PressureGrad)
    -0.15487679* (NormalStressGradV-PressureGrad)  0.44327040*(NormalStressGradV-PressureGrad)]
    RotationMatrix = 
    [cosd(RotationAngle) -sind(RotationAngle)
    sind(RotationAngle) cosd(RotationAngle)]
    RotatedStress = RotationMatrix * StressGradTensor * RotationMatrix'
    
    Friction_Main = -StressGradTensor[1,2]/StressGradTensor[1,1]
    Friction_Nucleation = Friction_Main + 0.01
    Friction_Branch= -RotatedStress[1,2]/RotatedStress[1,1]
    NormalStressGradMain = StressGradTensor[1,1]
    NormalStressGradBranch = RotatedStress[1,1]

    Theta_I = 1e10
    V_I = 1e-12
    V_0 = 1e-9
    Mu0 = Friction_Main - RSF_b * log(Theta_I * V_0 /RSF_Dc) - RSF_a * log(V_I/V_0)

    V_I_Nucleation = 1e-6
    Theta_I_Nucleation = RSF_Dc / V_0 * exp((Friction_Nucleation - Mu0 - RSF_a * log(V_I_Nucleation/V_0))/RSF_b)
    Theta_I_Branch = RSF_Dc / V_0 * exp((Friction_Branch - Mu0 - RSF_a * log(V_I/V_0))/RSF_b)





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
    Fault[1,1]=-8000.0
    Fault[1,2]=-0.0;
    Fault[1,3]=15000/2;
    Fault[1,4]=3e3;
    Fault[1,5]=3e3;
    Fault[1,6]=Fault1StrikeAngle;
    Fault[1,7]=Fault1DipAngle;

    ### Fault 2     
    Fault[2,1]=6e3
    Fault[2,2]=-0.0;
    Fault[2,3]=15000/2;
    Fault[2,4]=12000.;
    Fault[2,5]=Fault[1,3]*2;
    Fault[2,6]=Fault1StrikeAngle;
    Fault[2,7]=Fault1DipAngle;
    
    ### Fault 3     
    Fault[3,1]=12000 * cosd(30)/2
    Fault[3,2]=-12000 * sind(30)/2;
    Fault[3,3]=15000/2;
    Fault[3,4]=12000.;
    Fault[3,5]=Fault[2,3]*2;
    Fault[3,6]=Fault2StrikeAngle;
    Fault[3,7]=Fault2DipAngle;

    ### Fault 4.     
    Fault[4,1]=-8000.0
    Fault[4,2]=-0.0;
    Fault[4,3]=3e3;
    Fault[4,4]=3e3;
    Fault[4,5]=6e3;
    Fault[4,6]=Fault1StrikeAngle;
    Fault[4,7]=Fault1DipAngle;

    ### Fault 5.     
    Fault[5,1]=-8000.0
    Fault[5,2]=-0.0;
    Fault[5,3]=12e3;
    Fault[5,4]=3e3;
    Fault[5,5]=6e3;
    Fault[5,6]=Fault1StrikeAngle;
    Fault[5,7]=Fault1DipAngle;

    ### Fault 6.     
    Fault[6,1]=-3250.0
    Fault[6,2]=-0.0;
    Fault[6,3]=7.5e3;
    Fault[6,4]=6.5e3;
    Fault[6,5]=15e3;
    Fault[6,6]=Fault1StrikeAngle;
    Fault[6,7]=Fault1DipAngle;

    ### Fault 7.     
    Fault[7,1]=-3250.0 - 9.5e3
    Fault[7,2]=-0.0;
    Fault[7,3]=7.5e3;
    Fault[7,4]=6.5e3;
    Fault[7,5]=15e3;
    Fault[7,6]=Fault1StrikeAngle;
    Fault[7,7]=Fault1DipAngle;


    ### Common Values 
    Fault[:,9] .= RSF_a;
    Fault[:,10] .= RSF_b;
    Fault[:,11] .= RSF_Dc;
    Fault[:,12] .= Theta_I;
    Fault[:,13] .= V_I;

    ### Friction and NormalSTress
    
    ############################# Write Bulk Input #################################
    ######++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++######
    ## Bulk File Order 
    ##  1.Ctr_X     2.Ctr_Y 3.Ctr_Z 4.St_L	    5.Dip_L	    6.StAng	    7.DipAng	8.LR/RN (-1: LL / +1 RL or -1:Reverse / +1 Normal)
    ##  9.a         10.b	11.Dc	12.Theta_i	13. V_i     14. Friction_i 15.NormalStress at surface [Pa]  
    ##  16. NoarmalStress Gradient [Pa] 17. V_Const     18. Minimum Segment Length

    Fault[:,8].= 1
    Fault[:,14] .= Friction_Main
    Fault[:,15] .= NormalStressI #.* NormalStressParameter;
    Fault[:,16] .= NormalStressGradMain
    Fault[:,17] .=0.0
    Fault[:,18] .=MaxDiscritLength
    
    Fault[1,12] = Theta_I_Nucleation;
    Fault[1,13] = V_I_Nucleation;
    Fault[1,14] = Friction_Nucleation

    Fault[3,12] = Theta_I_Branch;
    Fault[3,14] = Friction_Branch
    Fault[3,16] = NormalStressGradBranch



    
    
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
        write(io, "Ctr_X\tCtr_Y\tCtr_Z\tSt_L\tDip_L\tStAng\tDipAng\tLR/RN\ta\tb\tDc\tTheta_i\tV_i\tFric_i\tSig0\tSigGrad\tV_Const\tMaxLeng\n")
        writedlm(io, InputFile)
    end;

    println("Saved File Name: ",OutputFileName)
end


BuildBulk()