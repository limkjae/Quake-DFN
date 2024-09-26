
using DelimitedFiles
using PyPlot
using PyCall
pygui(true)
include("../Functions_BuildInputFile.jl")
include("../Results/Functions_Plot.jl")

MaxStressOrientation = 135. # between 0-180 degree
StressRatioMaxOverMin = 0.5 # Sig3/Sig1
MaxSigSurface = 2e6 # Stress at surface
MaxSigGrad = 6000.0 # Stress gradient along max stress orientation [Pa/m]

Fault_Theta_i = 1e10 # if zero, automatically calculated based on the velocity and Mu0
Fault_V_i = 0.0 # if zero, automatically calculated based on the theta and Mu0
Friction_0 =  0.30 # if zero, automatically calculated based on the velocity and Theta

V0=1e-9
MinFrictionAllowed = 0.1


function adjustParameters(MaxStressOrientation, StressRatioMaxOverMin, MaxSigSurface, MaxSigGrad,  
                        Fault_Theta_i, Fault_V_i, Friction_0, V0, MinFrictionAllowed)

    InputBulkFileName="Input_BulkFaultGeometry.txt"

    Input_Bulk=readdlm(InputBulkFileName)
    Input_BulktoAdjust=Input_Bulk[4:end,:]

    ### Friction and NormalSTress
    for i=1:length(Input_BulktoAdjust[:,1])
        if Input_BulktoAdjust[i,17] == 0
            NormalStressParameter = (1+StressRatioMaxOverMin)/2 + (1-StressRatioMaxOverMin)/2 * cosd(2 * (Input_BulktoAdjust[i,6] - 90.0 - MaxStressOrientation))
            ShearStressParameter = -(1-StressRatioMaxOverMin)/2 * sind.(2 * (Input_BulktoAdjust[i,6] - 90.0 - MaxStressOrientation))
            Input_BulktoAdjust[i,8] = sign(ShearStressParameter)
            Input_BulktoAdjust[i,14] = abs(ShearStressParameter / NormalStressParameter)
            if Input_BulktoAdjust[i,14] < MinFrictionAllowed
                Input_BulktoAdjust[i,14] = MinFrictionAllowed
            end
            Input_BulktoAdjust[i,15] = MaxSigSurface * NormalStressParameter;
            Input_BulktoAdjust[i,16] = MaxSigGrad * NormalStressParameter
        end
    end

    

    if iszero(Fault_V_i)
        println("Initial Velocity adjusted")
        Input_BulktoAdjust[:,12] .= Fault_Theta_i
        Input_BulktoAdjust[:,13] = V0 .* exp.( (Input_BulktoAdjust[:,14] .- Friction_0 .- Input_BulktoAdjust[:,10] .* log.(Fault_Theta_i .* V0./Input_BulktoAdjust[:,11] )) ./ Input_BulktoAdjust[:,9] )
    end
    
    
    if iszero(Fault_Theta_i)
        println("Initial Theta adjusted")
        Input_BulktoAdjust[:,13] .= Fault_V_i
        Input_BulktoAdjust[:,12] = Input_BulktoAdjust[:,11] ./ V0 .* exp.( (Input_BulktoAdjust[:,14] .- Friction_0 .- Input_BulktoAdjust[:,9] .* log.(Input_BulktoAdjust[:,13] ./ V0)) ./ Input_BulktoAdjust[:,10])

    end
     
    if iszero(Friction_0)
        println("Theta and Velocity directly assigned")
        Input_BulktoAdjust[:,12] .= Fault_Theta_i
        Input_BulktoAdjust[:,13] .= Fault_V_i
    end

    ############################## Save File #############################
    
    open(InputBulkFileName, "w") do io
        write(io,"SwitchSS/RN\tShearMod\tPoissonRatio\tR_Density\tCrit_TooClose\tTooCloseNormal_Multiplier\tMinimum_NS\n")
        writedlm(io,[Input_Bulk[2,1]   Input_Bulk[2,2]     Input_Bulk[2,3]      Input_Bulk[2,4]   Input_Bulk[2,5]      Input_Bulk[2,6]  Input_Bulk[2,7] ])
        write(io, "Ctr_X\tCtr_Y\tCtr_Z\tSt_L\tDip_L\tStAng\tDipAng\tLR\ta\tb\tDc\tTheta_i\tV_i\tFric_i\tSig0\tSigGrad\tV_Const\tMaxLeng\n")
        writedlm(io, Input_BulktoAdjust)



#####################################  Plot  ##########################################
############################### Figure Configuration ##################################

PlotRotation=[60,-50]
Transparent = 1 # 1 for transparent fault plot
Edge = 1 # 0 for no element boudary 
MinMax_Axis=0
LoadingFaultCount=0 # 1 to plot constant velocity faults. 

PlotInput=Input_BulktoAdjust[:,8]; ColorMinMax=0 
# PlotInput= log10.(Input_BulktoAdjust[:,13]); ColorMinMax=0
# PlotInput=Input_BulktoAdjust[:,8] .* log10.(Input_BulktoAdjust[:,13]); ColorMinMax=0
# PlotInput= Result_NormalStress[PlotStep,:] -  Fault_NormalStress; ColorMinMax=[-1e6,1e6]
# PlotInput=ResultDisp[PlotStep,:]; ColorMinMax=0 

##### Input Code
##### 1.Ctr_X     2.Ctr_Y 3.Ctr_Z 4.St_L	    5.Dip_L	    6.StAng	    7.DipAng	8.LR/RN
#####  9.a         10.b	    11.Dc	12.Theta_i	13. V_i     14. Friction_i 15.NormalStress at surface [Pa]  
#####  16. NoarmalStress Gradient [Pa] 17. V_Const     18. Minimum Segment Length


figure(1)
clf()
MaxVaule, MinValue = FaultPlot_3D_Color_General(Input_BulktoAdjust[:,1:3],
Input_BulktoAdjust[:,4], Input_BulktoAdjust[:,5], Input_BulktoAdjust[:,6], Input_BulktoAdjust[:,7], Input_BulktoAdjust[:,8], PlotInput, 
    PlotRotation, MinMax_Axis, ColorMinMax, Transparent, Edge, LoadingFaultCount)

    # MaxVaule, MinValue = FaultPlot_3D_Color_General(FaultCenter,FaultLengthStrike, FaultLengthDip,
    # FaultStrikeAngle, FaultDipAngle, FaultLLRR, PlotInput, 
    # PlotRotation, MinMax_Axis, ColorMinMax, Transparent, Edge, LoadingFaultCount)

    ax = subplot(projection="3d")
    xlabel("x")
    ylabel("y")
plotforcbar=  scatter([1,1],[1,1],0.1, [MinValue,MaxVaule], cmap="jet")
colorbar(plotforcbar, pad=0.15)
figure(1).canvas.draw()


##############+++++++++++++++++++++++++++++++++++++++++++++++++++++++++################
#######################################################################################




    end
end

adjustParameters(MaxStressOrientation, StressRatioMaxOverMin, MaxSigSurface, MaxSigGrad,  
                Fault_Theta_i, Fault_V_i, Friction_0, V0, MinFrictionAllowed)

