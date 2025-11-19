
using DelimitedFiles
using PyPlot
using PyCall
using LinearAlgebra
using Statistics
pygui(true)
include("../scripts/Functions_BuildInputFile.jl")
include("../scripts/Functions_Plot.jl")
include("../scripts/Functions_StressApplication.jl")



InputBulkFileName="Input_BulkFaultGeometry.txt"

function ChangeBulk()

    ##############################################################################################
    ######################################## Inputs ##############################################
    ###### build Principal Stress. Compression Positive. Only Ratio Matters! ########
    PrincipalStressRatioX = 0.5
    PrincipalStressRatioY = 1.0
    PrincipalStressRatioZ = 0.3
    StressRotationStrike = 150 # degree
    StressRotationDip = -40   # degree

    MaximumTargetVelocity = 1e-11 # if this has value, the maximum velocity is set to this value. And Mu0 will be adjusted accordingly.
    ConstantMu0 = 0.0
    ConstantTheta = 1e10 # if not zero, initial theta will be revised to this uniformly
    Fault_a = 0 # if not zero, RSF "a" value will be revised to this uniformly
    Fault_b = 0 # if not zero, RSF "b" value will be revised to this uniformly
    Fault_Dc = 1e-3 # if not zero, RSF "Dc" value will be revised to this uniformly

    MinFrictionAllowed = 0.05
    ShearModulus = 30e9
    PoissonRatio = 0.25
    Rock_Density = 2700.0
    MinimumNormalStressAllowed = 1e6
    StressOnSurface_Sig1Orientation = 20e6 # pascal
    StressGredient_Sig1Orientation = 0 # pascal/m

    FaultSegmentLength = 0 # if 0, segment length will be unchanged (Only For Rectangular Grid)   
    LoadingFaultInvert = 1 # if 1, loading fault sense of slip become inverted



    ############################  FigureConfiguration  ##############################
    PlotProperty = 1 # 1: plot traction, normal, shear vector, 0: no
    PlotSenseofSlip = 1 # 1: plot loading fault, normal, shear vector, 0: no
    # PlotSenseofSlipArrow = 1 # 1: plot arrow for sense of slip, 0: no
    PlotLoadingFault = 0


    ##### Which Property? 
    ##### 1: Rake	2:a	3:b	4:Dc 5:	Theta_i	6:V_i	7: Fric_i	8: Sig0	9:SigGrad	10:V_Const	########
    WhattoPlot = 7

    PlotRotation=[63,-120]
    Transparent = 1 # 1 for transparent fault plot
    Edge = 1 # 0 for no element boudary 
    MinMax_Axis=0

    StressVectorLocation = 0 # Autometically Adjusted when 0 
    PrinpalStressLength = 100 # Autometically Adjusted when 0 
    V_p = 1e-5 # When target velocity is set this will be used for peak friction plot
    V_r = 1e-2 # When target velocity is set this will be used for residual friction plot
    #################################################################################
    ##############################################################################################
    ##############################################################################################



    ####################### Principal stress process #########################

    MaxStressRatio =  maximum([PrincipalStressRatioX,PrincipalStressRatioY, PrincipalStressRatioZ])
    StressRotationStrike = StressRotationStrike + 0.001
    StressRotationDip = StressRotationDip + 0.001
    PrincipalStressRatioX = PrincipalStressRatioX / MaxStressRatio
    PrincipalStressRatioY = PrincipalStressRatioY / MaxStressRatio
    PrincipalStressRatioZ = PrincipalStressRatioZ / MaxStressRatio
    Sig1 = maximum([PrincipalStressRatioX, PrincipalStressRatioY, PrincipalStressRatioZ])
    Sig3 = minimum([PrincipalStressRatioX, PrincipalStressRatioY, PrincipalStressRatioZ])
        
    #################### Calculate Stress for XYZ Frame #####################
    # Stres Negative for Compression
    PrincipalStressRatio = [-PrincipalStressRatioX 0 0 
                            0 -PrincipalStressRatioY 0
                            0 0 -PrincipalStressRatioZ]
    # StressRotationStrike = -StressRotationStrike
    # Rotate the Stress along Maximum Stress Angle (XYZ corrdinate)
    RotationMat_Strike =
    [cosd(StressRotationStrike) -sind(StressRotationStrike)  0
    sind(StressRotationStrike) cosd(StressRotationStrike) 0
    0  0  1]

    RotationMat_Dip =
    [1 0 0
    0 cosd(StressRotationDip) -sind(StressRotationDip) 
    0 sind(StressRotationDip) cosd(StressRotationDip)]

    RotationMat =  RotationMat_Strike * RotationMat_Dip
    StressRatioXYZ = RotationMat * PrincipalStressRatio * RotationMat'  
    ########^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^########


    ######################### Read Bulk Fault Geometry ######################
    Input_Bulk=readdlm(InputBulkFileName)
    Input_Bulk_Header = Input_Bulk[2,:]
    Switch_StrikeSlip_or_ReverseNormal = Input_Bulk[2,1] 
    Input_Bulk_Header[1] = 0.0 # SwitchSS/RN
    if ShearModulus !=0; Input_Bulk_Header[2] = ShearModulus; end
    if PoissonRatio !=0; Input_Bulk_Header[3] = PoissonRatio; end
    if Rock_Density !=0; Input_Bulk_Header[4] = Rock_Density; end
    if MinimumNormalStressAllowed !=0; Input_Bulk_Header[7] = MinimumNormalStressAllowed; end
    Input_Bulk_Header = filter(!isempty, Input_Bulk_Header)

    Input_Bulk=Input_Bulk[4:end,:]
    BulkFaultCount = length(Input_Bulk[:,1])
    if size(Input_Bulk, 2) == 18
        println("Rectangle")
        RorT = "R"
        Input_Bulk = LRtoRake(Switch_StrikeSlip_or_ReverseNormal, Input_Bulk)
        Input_Bulk = Input_Bulk[sortperm(Input_Bulk[:, 17], rev=false), :] # move the loading faults to the end
        LoadingFaultCount = sum(x->x>0, Input_Bulk[:,17])
    elseif  size(Input_Bulk, 2) == 20
        println("Triangle")
        RorT = "T"
        Input_Bulk = Input_Bulk[sortperm(Input_Bulk[:, 19], rev=false), :] # move the loading faults to the end
        LoadingFaultCount = sum(x->x>0, Input_Bulk[:,19])
    else 
        error("Input Bulk Fault Geometry file should have 18 or 20 columns")
    end

    ########################################################################



    ##################### Define parameters for figures #########################
    if PlotLoadingFault == 0 
        LoadingFaultCountPlot = LoadingFaultCount
    else 
        LoadingFaultCountPlot = 0
    end

    if StressVectorLocation == 0
        StressVectorLocation = [mean(Input_Bulk[1:end - LoadingFaultCountPlot,1]), mean(Input_Bulk[1:end - LoadingFaultCountPlot,2]), 100]
    end 
    if PrinpalStressLength == 0
        PrinpalStressLength = maximum([maximum(Input_Bulk[1:end-LoadingFaultCountPlot,1]) - minimum(Input_Bulk[1:end-LoadingFaultCountPlot,1]), 
                                      maximum(Input_Bulk[1:end-LoadingFaultCountPlot,2]) - minimum(Input_Bulk[1:end-LoadingFaultCountPlot,2])])/2
    end





    ############################# Calcuate Stress and Frictions in Each Fault ##########################

    if RorT == "R"
        Input_Bulk = CalculateFrictionStress_R(BulkFaultCount, Input_Bulk, MinFrictionAllowed,StressRatioXYZ,LoadingFaultInvert,
                        StressOnSurface_Sig1Orientation, StressGredient_Sig1Orientation, MinimumNormalStressAllowed)
        if Fault_a != 0.0; Input_Bulk[:,9] .= Fault_a; end
        if Fault_b != 0.0; Input_Bulk[:,10] .= Fault_b; end
        if Fault_Dc != 0.0; Input_Bulk[:,11] .= Fault_Dc; end
        if ConstantTheta != 0.0; Input_Bulk[:,12] .= ConstantTheta; end
    elseif RorT == "T"
        Input_Bulk, UnitVector_Normal, UnitVector_Slip, UnitVector_DipSlip, UnitVector_StrikeSlip = 
                        CalculateFrictionStress_T(BulkFaultCount, Input_Bulk, MinFrictionAllowed,StressRatioXYZ, 
                        StressOnSurface_Sig1Orientation, StressGredient_Sig1Orientation, MinimumNormalStressAllowed, LoadingFaultCount,LoadingFaultInvert) 

        if Fault_a != 0.0; Input_Bulk[:,11] .= Fault_a; end
        if Fault_b != 0.0; Input_Bulk[:,12] .= Fault_b; end
        if Fault_Dc != 0.0; Input_Bulk[:,13] .= Fault_Dc; end
        if ConstantTheta != 0.0; Input_Bulk[:,14] .= ConstantTheta; end
    end
    if maximum(Input_Bulk[:,16]) > 1
        println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        println("!!!!!!! Maximum Friction is larger than 1. Check the input parameters. Maximum Friction: ", maximum(Input_Bulk[:,16]))
        println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    end





    ####################################### Adjust Velocity #############################################


    if RorT == "R"
    Input_Bulk, Mu0, PeakFriction, ResidualFriction = 
        FindMu0_AdjV_R(Input_Bulk, MaximumTargetVelocity, V_p, V_r, ConstantMu0)
        if FaultSegmentLength > 0
            Input_Bulk[:,18]  .= FaultSegmentLength

        end
    elseif RorT == "T"
    Input_Bulk, Mu0, PeakFriction, ResidualFriction = 
        FindMu0_AdjV_T(Input_Bulk, MaximumTargetVelocity, V_p, V_r, ConstantMu0)
    end



    #################### Remove Undefined Rake Case ######################
    
    if RorT == "R"
        Input_BulkFiltered = zeros(1,18) 
        RakeIndex = 8
    elseif RorT == "T"
        Input_BulkFiltered = zeros(1,20) 
        RakeIndex = 10
    end

    SurvivedFaults = 0
    for i=1:BulkFaultCount
        SurvivedFaults =+ 1
        if isnan(Input_Bulk[i,RakeIndex])
            println("Slip sence of Bulk Element ",i, " cannot be defined with the given stress field. The element will be removed")
            println("This problem can be alleviated by adding anisotropy or small angles to the stress field applied")
        else
            Input_BulkFiltered = [Input_BulkFiltered; Input_Bulk[i,:]']
        end
    end
    Input_BulkFiltered = Input_BulkFiltered[2:end,:]


    ############################## Save File #############################

    open(InputBulkFileName, "w") do io
        write(io,"SwitchSS/RN\tShearMod\tPoissonRatio\tR_Density\tCrit_TooClose\tTooCloseNormal_Multiplier\tMinimum_NS\n")
        writedlm(io, Input_Bulk_Header')
        write(io, "Ctr_X\tCtr_Y\tCtr_Z\tSt_L\tDip_L\tStAng\tDipAng\tRake\ta\tb\tDc\tTheta_i\tV_i\tFric_i\tSig0\tSigGrad\tV_Const\tMaxLeng\n")
        writedlm(io, Input_BulkFiltered)
    end


    





        ######################################### Plotting #############################################
    if PlotProperty == 1
        ######################################## Surface Plot ##########################################
        if RorT == "R"
            figure(1)
            clf()
            if WhattoPlot == 5 || WhattoPlot == 6  
                PlotInput = log10.(Input_Bulk[1:end-LoadingFaultCountPlot,WhattoPlot+7]); ColorMinMax = 0
            else
                PlotInput = Input_Bulk[1:end-LoadingFaultCountPlot,WhattoPlot+7]; ColorMinMax = 0
            end
            MaxValue, MinValue = FaultPlot_3D_Color_General(Input_Bulk[1:end-LoadingFaultCountPlot,1:3],
                Input_Bulk[1:end-LoadingFaultCountPlot,4], Input_Bulk[1:end-LoadingFaultCountPlot,5], Input_Bulk[1:end-LoadingFaultCountPlot,6], 
                Input_Bulk[1:end-LoadingFaultCountPlot,7], Input_Bulk[1:end-LoadingFaultCountPlot,8], PlotInput, 
                PlotRotation, MinMax_Axis, ColorMinMax, Transparent, Edge, 0)
                ax = subplot(projection="3d")
                xlabel("x")
                ylabel("y")
                plotforcbar=  scatter([1,1],[1,1],0.1, [MinValue,MaxValue], cmap="jet")
                cbar  = colorbar(plotforcbar, pad=0.15)
                figure(1).canvas.draw()
        elseif RorT == "T"
            if WhattoPlot == 5 || WhattoPlot == 6  
                PlotInput = log10.(Input_Bulk[1:end-LoadingFaultCountPlot,WhattoPlot+9]); ColorMinMax = 0
            else
                PlotInput = Input_Bulk[1:end-LoadingFaultCountPlot,WhattoPlot+9]; ColorMinMax = 0
            end
            # InputProperty = log10.(Input_Bulk[1:end-LoadingFaultCountPlot,15]); ColorMinMax = 0
            MaxValue=maximum(PlotInput) * 1.001
            MinValue=minimum(PlotInput) * 0.999

            ArrowLength = 500
            figure(1)
            clf()
            art3d = PyObject(PyPlot.art3D)
            ax = subplot(projection="3d")
            for ElemIdx = 1:BulkFaultCount - LoadingFaultCountPlot
                cm = get_cmap(:jet)
                PlotValue=(PlotInput[ElemIdx]-MinValue)/(MaxValue-MinValue)
                # println(PlotValue)
                face_color = [cm(PlotValue)[1], cm(PlotValue)[2],cm(PlotValue)[3], 0.5]

                verts = ((Input_Bulk[ElemIdx,1:3],Input_Bulk[ElemIdx,4:6],Input_Bulk[ElemIdx,7:9]), )
                p3c = PyObject(art3d.Poly3DCollection(verts))
                pycall(ax.add_collection3d, PyAny, p3c)

                edge_color = [0.2, 0.2, 0.2, 0.5]

                pycall(p3c.set_facecolor, PyAny, face_color)
                pycall(p3c.set_edgecolor, PyAny, edge_color)
                ax.view_init(45, -30)

            end
                plotforcbar=  scatter([1,1],[1,1],0.1, [MinValue,MaxValue], cmap="jet")
                
                # plotforcbar=  scatter([1,1,1],[1,1,1],[1,1,1],0.1, [MinValue,MaxValue], cmap="jet")

                colorbar(plotforcbar, pad=0.15)
            xlabel("x")
            ylabel("y")
            ax.set_aspect("equal")

        end

        ##################################### Principal Stresses #######################################
            Linewidth = 2
            Arrow_length_ratio = 0.01

            UnlotatedVectorX =[PrincipalStressRatioX, 0.0, 0.0]
            UnlotatedVectorY =[0.0, PrincipalStressRatioY, 0.0]
            UnlotatedVectorZ =[0.0, 0.0, PrincipalStressRatioZ]

            RotationMat_Strike=
            [cosd(StressRotationStrike) -sind(StressRotationStrike)  0
            sind(StressRotationStrike) cosd(StressRotationStrike) 0
            0  0  1];

            RotationMat_Dip=
            [1 0 0
            0 cosd(StressRotationDip) -sind(StressRotationDip)
            0 sind(StressRotationDip) cosd(StressRotationDip)]

            RotatedStessX = RotationMat_Strike * RotationMat_Dip  * UnlotatedVectorX * PrinpalStressLength
            RotatedStessY = RotationMat_Strike * RotationMat_Dip  * UnlotatedVectorY * PrinpalStressLength
            RotatedStessZ = RotationMat_Strike * RotationMat_Dip  * UnlotatedVectorZ * PrinpalStressLength

            ax.quiver(StressVectorLocation[1] + RotatedStessX[1], StressVectorLocation[2] + RotatedStessX[2], StressVectorLocation[3]+ RotatedStessX[3], 
            -RotatedStessX[1] , -RotatedStessX[2] , -RotatedStessX[3] ,
            color="k",arrow_length_ratio=Arrow_length_ratio, linewidth =Linewidth)

            ax.quiver(StressVectorLocation[1] + RotatedStessY[1], StressVectorLocation[2] + RotatedStessY[2], StressVectorLocation[3]+ RotatedStessY[3], 
            -RotatedStessY[1] , -RotatedStessY[2] , -RotatedStessY[3] ,
            color="k",arrow_length_ratio=Arrow_length_ratio, linewidth =Linewidth)

            ax.quiver(StressVectorLocation[1] + RotatedStessZ[1], StressVectorLocation[2] + RotatedStessZ[2], StressVectorLocation[3]+ RotatedStessZ[3], 
            -RotatedStessZ[1] , -RotatedStessZ[2] , -RotatedStessZ[3] ,
            color="k",arrow_length_ratio=Arrow_length_ratio, linewidth =Linewidth)

            ax.quiver(StressVectorLocation[1] - RotatedStessX[1], StressVectorLocation[2] - RotatedStessX[2], StressVectorLocation[3] - RotatedStessX[3], 
            RotatedStessX[1] , RotatedStessX[2] , RotatedStessX[3] ,
            color="k",arrow_length_ratio=Arrow_length_ratio, linewidth =Linewidth)

            ax.quiver(StressVectorLocation[1] - RotatedStessY[1], StressVectorLocation[2] - RotatedStessY[2], StressVectorLocation[3] - RotatedStessY[3], 
            RotatedStessY[1] , RotatedStessY[2] , RotatedStessY[3] ,
            color="k",arrow_length_ratio=Arrow_length_ratio, linewidth =Linewidth)

            ax.quiver(StressVectorLocation[1] - RotatedStessZ[1], StressVectorLocation[2] - RotatedStessZ[2], StressVectorLocation[3] - RotatedStessZ[3], 
            RotatedStessZ[1] , RotatedStessZ[2] , RotatedStessZ[3] ,
            color="k",arrow_length_ratio=Arrow_length_ratio, linewidth =Linewidth)

            ax.set_aspect("equal")
            
    end
    if PlotSenseofSlip == 1

        ##################################### Sense of Slip #######################################

        if RorT == "R"
            figure(2)
            clf()

            ax = PlotBulk_SenseOfSlip(0.0, Input_Bulk[1:end-LoadingFaultCountPlot,:], PlotRotation, Transparent, Edge, MinMax_Axis)


        elseif RorT == "T"
            RakeAngle_NRAdjusted = copy(Input_Bulk[1:end-LoadingFaultCountPlot,10])
            for BulkIndex = 1: BulkFaultCount - LoadingFaultCountPlot
                if UnitVector_Normal[BulkIndex,3] > 0
                    RakeAngle_NRAdjusted[BulkIndex] = 360 - RakeAngle_NRAdjusted[BulkIndex]

                end
            end
            PlotInput=RakeAngle_NRAdjusted
            MaxValue = 360.0
            MinValue = 0.0
            figure(2)
            clf()
            art3d = PyObject(PyPlot.art3D)
            ax = subplot(projection="3d")
            for ElemIdx = 1:BulkFaultCount - LoadingFaultCountPlot
                cm = get_cmap(:hsv)
                PlotValue=(PlotInput[ElemIdx]-MinValue)/(MaxValue-MinValue)

                face_color = [cm(PlotValue)[1], cm(PlotValue)[2],cm(PlotValue)[3], 0.5]

                verts = ((Input_Bulk[ElemIdx,1:3],Input_Bulk[ElemIdx,4:6],Input_Bulk[ElemIdx,7:9]), )
                p3c = PyObject(art3d.Poly3DCollection(verts))
                pycall(ax.add_collection3d, PyAny, p3c)

                edge_color = [0.2, 0.2, 0.2, 0.5]

                pycall(p3c.set_facecolor, PyAny, face_color)
                pycall(p3c.set_edgecolor, PyAny, edge_color)
                ax.view_init(PlotRotation[1], PlotRotation[2])

            end
            plotforcbar=  scatter([1,1],[1,1],0.1, [MinValue,MaxValue], cmap="hsv")
            cbar  = colorbar(plotforcbar, pad=0.15)
            cbar.set_ticks([0, 90, 180,270,360])
            cbar.set_ticklabels(["Left Lateral", "Reverse", "Right Lateral", "Normal", "Left Lateral"])
            xlabel("x")
            ylabel("y")
            ax.set_aspect("equal")
        end
            ##################################### Principal Stresses #######################################
            Linewidth = 2
            Arrow_length_ratio = 0.01

            UnlotatedVectorX =[PrincipalStressRatioX, 0.0, 0.0]
            UnlotatedVectorY =[0.0, PrincipalStressRatioY, 0.0]
            UnlotatedVectorZ =[0.0, 0.0, PrincipalStressRatioZ]

            RotationMat_Strike=
            [cosd(StressRotationStrike) -sind(StressRotationStrike)  0
            sind(StressRotationStrike) cosd(StressRotationStrike) 0
            0  0  1];

            RotationMat_Dip=
            [1 0 0
            0 cosd(StressRotationDip) -sind(StressRotationDip)
            0 sind(StressRotationDip) cosd(StressRotationDip)]

            RotatedStessX = RotationMat_Strike * RotationMat_Dip  * UnlotatedVectorX * PrinpalStressLength
            RotatedStessY = RotationMat_Strike * RotationMat_Dip  * UnlotatedVectorY * PrinpalStressLength
            RotatedStessZ = RotationMat_Strike * RotationMat_Dip  * UnlotatedVectorZ * PrinpalStressLength

            ax.quiver(StressVectorLocation[1] + RotatedStessX[1], StressVectorLocation[2] + RotatedStessX[2], StressVectorLocation[3]+ RotatedStessX[3], 
            -RotatedStessX[1] , -RotatedStessX[2] , -RotatedStessX[3] ,
            color="k",arrow_length_ratio=Arrow_length_ratio, linewidth =Linewidth)

            ax.quiver(StressVectorLocation[1] + RotatedStessY[1], StressVectorLocation[2] + RotatedStessY[2], StressVectorLocation[3]+ RotatedStessY[3], 
            -RotatedStessY[1] , -RotatedStessY[2] , -RotatedStessY[3] ,
            color="k",arrow_length_ratio=Arrow_length_ratio, linewidth =Linewidth)

            ax.quiver(StressVectorLocation[1] + RotatedStessZ[1], StressVectorLocation[2] + RotatedStessZ[2], StressVectorLocation[3]+ RotatedStessZ[3], 
            -RotatedStessZ[1] , -RotatedStessZ[2] , -RotatedStessZ[3] ,
            color="k",arrow_length_ratio=Arrow_length_ratio, linewidth =Linewidth)

            ax.quiver(StressVectorLocation[1] - RotatedStessX[1], StressVectorLocation[2] - RotatedStessX[2], StressVectorLocation[3] - RotatedStessX[3], 
            RotatedStessX[1] , RotatedStessX[2] , RotatedStessX[3] ,
            color="k",arrow_length_ratio=Arrow_length_ratio, linewidth =Linewidth)

            ax.quiver(StressVectorLocation[1] - RotatedStessY[1], StressVectorLocation[2] - RotatedStessY[2], StressVectorLocation[3] - RotatedStessY[3], 
            RotatedStessY[1] , RotatedStessY[2] , RotatedStessY[3] ,
            color="k",arrow_length_ratio=Arrow_length_ratio, linewidth =Linewidth)

            ax.quiver(StressVectorLocation[1] - RotatedStessZ[1], StressVectorLocation[2] - RotatedStessZ[2], StressVectorLocation[3] - RotatedStessZ[3], 
            RotatedStessZ[1] , RotatedStessZ[2] , RotatedStessZ[3] ,
            color="k",arrow_length_ratio=Arrow_length_ratio, linewidth =Linewidth)

            ax.set_aspect("equal")
            


            ################# Draw Sense of Slip ################
            RakeAngle = copy(Input_Bulk[1:end-LoadingFaultCountPlot,10])
            FaultCenter = (Input_Bulk[:,1:3] + Input_Bulk[:,4:6] + Input_Bulk[:,7:9])/3
            FaultSize = mean([norm.(eachrow(Input_Bulk[:,1:3] - Input_Bulk[:,4:6]))  norm.(eachrow(Input_Bulk[:,1:3] - Input_Bulk[:,7:9]))  norm.(eachrow(Input_Bulk[:,7:9] - Input_Bulk[:,4:6]))], dims = 2)

            LineLength = FaultSize /4
            UnrotatedSlipUnitVec = [cosd.(RakeAngle) sind.(RakeAngle) zeros(BulkFaultCount-LoadingFaultCountPlot)]
            UnrotatedGapVector = [0 0 1]
            RotatedSlipUnitVec = UnrotatedSlipUnitVec .* 0.0 
            VecStart = UnrotatedSlipUnitVec .* 0.0 
            VecEnd = UnrotatedSlipUnitVec .* 0.0 
            

        PlotSenseofSlipArrow =0 
        if PlotSenseofSlipArrow == 1 
            for BulkIndex = 1: BulkFaultCount - LoadingFaultCountPlot

            
            RotatedGap = UnitVector_Normal[BulkIndex,:] .* FaultSize[BulkIndex ] * 0.1

            VecStart[BulkIndex,:] = -UnitVector_Slip[BulkIndex,:] /2 .* LineLength[BulkIndex] .+ FaultCenter[BulkIndex,1:3] - RotatedGap
            VecEnd[BulkIndex,:] = UnitVector_Slip[BulkIndex,:]/2 .* LineLength[BulkIndex] .+ FaultCenter[BulkIndex,1:3] + RotatedGap
            end

            for i =1:BulkFaultCount - LoadingFaultCountPlot
                ax.quiver(VecStart[i,1], VecStart[i,2], VecStart[i,3], 
                    UnitVector_Slip[i,1] * LineLength[i], UnitVector_Slip[i,2] * LineLength[i], UnitVector_Slip[i,3] * LineLength[i],
                    color="k",arrow_length_ratio=0.2)
                ax.quiver(VecEnd[i,1], VecEnd[i,2], VecEnd[i,3], 
                    -UnitVector_Slip[i,1] * LineLength[i], -UnitVector_Slip[i,2] * LineLength[i], -UnitVector_Slip[i,3] * LineLength[i],
                    color="k",arrow_length_ratio=0.2)
           
                ax.quiver(FaultCenter[i,1], FaultCenter[i,2], FaultCenter[i,3], 
                    UnitVector_Normal[i,1] * LineLength[i], UnitVector_Normal[i,2] * LineLength[i], UnitVector_Normal[i,3] * LineLength[i],
                    color="r",arrow_length_ratio=0.2)

                ax.quiver(FaultCenter[i,1], FaultCenter[i,2], FaultCenter[i,3], 
                    UnitVector_DipSlip[i,1] * LineLength[i], UnitVector_DipSlip[i,2] * LineLength[i], UnitVector_DipSlip[i,3] * LineLength[i],
                    color="b",arrow_length_ratio=0.2)
                    
                ax.quiver(FaultCenter[i,1], FaultCenter[i,2], FaultCenter[i,3], 
                    UnitVector_StrikeSlip[i,1] * LineLength[i], UnitVector_StrikeSlip[i,2] * LineLength[i], UnitVector_StrikeSlip[i,3] * LineLength[i],
                    color="g",arrow_length_ratio=0.2)

            end
        end
    end

end

ChangeBulk()

include("../Plot_BulkFaultGeometry.jl")