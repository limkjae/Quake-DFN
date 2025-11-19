

function lltokm(latlon1,latlon2)
    radius=6371.0;
    lon1=0*pi/180
    lon2=latlon1*pi/180
    lat1=0*pi/180
    lat2=latlon2*pi/180
    deltaLat=lat2-lat1
    deltaLon=lon2-lon1
    a=sin((deltaLat)/2)^2 + cos(lat1)*cos(lat2) * sin(deltaLon/2)^2
    #c=2*atan2(sqrt(a),sqrt(1-a))
    c=2*atan(sqrt(a)/sqrt(1-a))
    x=deltaLon*cos((lat1+lat2)/2)*radius
    y=deltaLat*radius
    return x, y
end


function FaultPlot_3D(FaultCenter,FaultLengthStrike, FaultLengthDip, FaultStrikeAngle, FaultDipAngle, FaultLLRR)
    
    cm = get_cmap(:jet)
    buffer=maximum(FaultLengthStrike)/3
    xmin=minimum(FaultCenter[:,1])-buffer
    xmax=maximum(FaultCenter[:,1])+buffer
    ymin=minimum(FaultCenter[:,2])-buffer
    ymax=maximum(FaultCenter[:,2])+buffer
    zmax=0
    zmin=-maximum(FaultCenter[:,3])-buffer
    xCenter=(xmin+xmax)/2
    yCenter=(ymin+ymax)/2
    zCenter=(zmin+zmax)/2
    
    maxlength=maximum([xmax-xmin,ymax-ymin, zmin-zmax])
    xmin=xCenter-maxlength/2
    xmax=xCenter+maxlength/2
    ymin=yCenter-maxlength/2
    ymax=yCenter+maxlength/2
    zmin=zmax-maxlength-buffer

    MaxValue=maximum(FaultLLRR)
    MinValue=minimum(FaultLLRR)

    for FaultIdx in eachindex(FaultLengthStrike)
        #print(FaultIdx)
        RotMatStrike=[cosd(FaultStrikeAngle[FaultIdx]) -sind(FaultStrikeAngle[FaultIdx]) 0
            sind(FaultStrikeAngle[FaultIdx]) cosd(FaultStrikeAngle[FaultIdx]) 0
            0  0 1]
        RotMatDip=[1 0  0
        0 cosd(FaultDipAngle[FaultIdx]) -sind(FaultDipAngle[FaultIdx])
        0 sind(FaultDipAngle[FaultIdx]) cosd(FaultDipAngle[FaultIdx])]
                
        p1=RotMatStrike*RotMatDip*[FaultLengthStrike[FaultIdx]/2;-FaultLengthDip[FaultIdx]/2;0] + [FaultCenter[FaultIdx,1]; FaultCenter[FaultIdx,2]; -FaultCenter[FaultIdx,3]];
        p2=RotMatStrike*RotMatDip*[-FaultLengthStrike[FaultIdx]/2;-FaultLengthDip[FaultIdx]/2;0] + [FaultCenter[FaultIdx,1]; FaultCenter[FaultIdx,2]; -FaultCenter[FaultIdx,3]];
        p3=RotMatStrike*RotMatDip*[-FaultLengthStrike[FaultIdx]/2;+FaultLengthDip[FaultIdx]/2;0]+ [FaultCenter[FaultIdx,1]; FaultCenter[FaultIdx,2]; -FaultCenter[FaultIdx,3]];
        p4=RotMatStrike*RotMatDip*[FaultLengthStrike[FaultIdx]/2;+FaultLengthDip[FaultIdx]/2;0]+ [FaultCenter[FaultIdx,1]; FaultCenter[FaultIdx,2]; -FaultCenter[FaultIdx,3]];
        
        fig = figure(3)
        #ion() # ioff() for turn off drawnow
        art3d = PyObject(PyPlot.art3D)

        verts2 = ([tuple(p1...); tuple(p2...); tuple(p3...); tuple(p4...)],)

        p3c = PyObject(art3d.Poly3DCollection(verts2, linewidths=1, alpha=0.5))
        ax = subplot(projection="3d")
        pycall(ax.add_collection3d, PyAny, p3c)
        xlim(xmin,xmax )
        ylim(ymin,ymax )
        zlim(zmin,zmax )
        PlotValue=(FaultLLRR[FaultIdx] - MinValue) / (MaxValue - MinValue)
        face_color = [cm(PlotValue)[1], cm(PlotValue)[2],cm(PlotValue)[3],1]
        
        # face_color = [0, 0, 1]
        edge_color = [0, 0, 0]
        pycall(p3c.set_facecolor, PyAny, face_color)
        pycall(p3c.set_edgecolor, PyAny, edge_color)
    end
   
end


function BulkToSegment(Input_Bulk)
    ## Bulk InputFile Order 
    ##  1.Ctr_X     2.Ctr_Y 3.Ctr_Z 4.St_L	    5.Dip_L	    6.StAng	    7.DipAng	8.LR
    ##  9.a         10.b	11.Dc	12.Theta_i	13. V_i     14. Friction_i 15.NormalStress  
    ##  16. NoarmalStress Gradient [Pa] 17. V_Const     18. Minimum Segment Length

    FaultCenter_Bulk=Input_Bulk[:,1:3]
    FaultLengthStrike_Bulk=Input_Bulk[:,4]
    FaultLengthDip_Bulk=Input_Bulk[:,5]
    FaultStrikeAngle_Bulk=Input_Bulk[:,6]
    FaultDipAngle_Bulk=Input_Bulk[:,7]
    Fault_MaximumSegmentLength=Input_Bulk[:,18]

    FaultSegmentCount=zeros(1,2)
    for i in eachindex(FaultLengthStrike_Bulk)    
        FaultSegmentCount=[FaultSegmentCount;[ceil(FaultLengthStrike_Bulk[i]/Fault_MaximumSegmentLength[i]),ceil(FaultLengthDip_Bulk[i]/Fault_MaximumSegmentLength[i])]']
    end
    FaultSegmentCount=FaultSegmentCount[2:end,:]
    TotalSegmentCount=sum(Int64,FaultSegmentCount[:,1].*FaultSegmentCount[:,2])
    println("Total element count is ", TotalSegmentCount)
    
    Input_Segment=zeros(TotalSegmentCount,19)

    ## Output Segment File Order 
    ##  1.Ctr_X     2.Ctr_Y 3.Ctr_Z 4.St_L	    5.Dip_L	    6.StAng	    7.DipAng	8.LR
    ##  9.a         10.b	11.Dc	12.Theta_i	13. V_i     14. Friction_i 15.NormalStress  
    ##  16. V_Const 17. Bulk Number     18. Bulk Strike Length      19. Bulk Dip Length

    FaultSegmentIdx=0;

    FaultSegmentIdx=0;
    FaultSegmentToBulk=0
    LengthSegmentCount=0
    WidthSegmentCount=0
    FaultStrikeAngle=0.0
    FaultDipAngle=0.0
    FaultLLRR=0.0
    FaultLengthStrike=0.0
    FaultLengthDip=0.0
    FaultCenter=zeros(1,3)

    for FaultBulkIdx in eachindex(FaultLengthStrike_Bulk)

        LengthSegmentCount=FaultSegmentCount[FaultBulkIdx,1];
        WidthSegmentCount=FaultSegmentCount[FaultBulkIdx,2];
        for LengthSegmentIdx=1:LengthSegmentCount
           for WidthSegmentIdx=1:WidthSegmentCount
                FaultSegmentIdx=FaultSegmentIdx+1;

                FaultLengthStrike=FaultLengthStrike_Bulk[FaultBulkIdx]/LengthSegmentCount;
                FaultLengthDip=FaultLengthDip_Bulk[FaultBulkIdx]/WidthSegmentCount;
                #println(FaultLengthStrike)



                FaultCenterB=FaultCenter_Bulk[FaultBulkIdx,:]
                FaultStrikeAngleB=FaultStrikeAngle_Bulk[FaultBulkIdx]
                FaultDipAngleB=FaultDipAngle_Bulk[FaultBulkIdx]
                FaultLengthStrikeB=FaultLengthStrike_Bulk[FaultBulkIdx]
                FaultLengthDipB=FaultLengthDip_Bulk[FaultBulkIdx]

                if LengthSegmentCount*WidthSegmentCount==1
                    FaultCenter=FaultCenterB;
                else
                    RotMatStrike=[cosd(FaultStrikeAngleB) -sind(FaultStrikeAngleB) 0
                    sind(FaultStrikeAngleB) cosd(FaultStrikeAngleB) 0
                    0  0 1]
                    RotMatDip=[1 0  0
                    0 cosd(FaultDipAngleB) -sind(FaultDipAngleB)
                    0 sind(FaultDipAngleB) cosd(FaultDipAngleB)]
                    
                    FaultCenter=RotMatStrike*RotMatDip* 
                    [-FaultLengthStrikeB/2 + FaultLengthStrikeB*(2*LengthSegmentIdx-1)/LengthSegmentCount/2;
                        -FaultLengthDipB/2 + FaultLengthDipB*(2*WidthSegmentIdx-1)/WidthSegmentCount/2 ; 0]+
                        [FaultCenterB[1]; FaultCenterB[2];  -FaultCenterB[3]]      
                    FaultCenter[3]=-FaultCenter[3]

                end

                Input_Segment[FaultSegmentIdx,1:3]=FaultCenter
                Input_Segment[FaultSegmentIdx,4]=FaultLengthStrike
                Input_Segment[FaultSegmentIdx,5]=FaultLengthDip
                Input_Segment[FaultSegmentIdx,6:14]=Input_Bulk[FaultBulkIdx,6:14]
                Input_Segment[FaultSegmentIdx,15]=Input_Bulk[FaultBulkIdx,15]+Input_Bulk[FaultBulkIdx,16]*FaultCenter[3]
                Input_Segment[FaultSegmentIdx,16]=Input_Bulk[FaultBulkIdx,17]
                Input_Segment[FaultSegmentIdx,17]=FaultBulkIdx
                Input_Segment[FaultSegmentIdx,18]=Input_Bulk[FaultBulkIdx,4]
                Input_Segment[FaultSegmentIdx,19]=Input_Bulk[FaultBulkIdx,5]


           end
        end
    end

    
    return Input_Segment

end


function LRtoRake(Switch_StrikeSlip_or_ReverseNormal, Input_Bulk)
    BulkFaultcount = size(Input_Bulk, 1)

    if Switch_StrikeSlip_or_ReverseNormal == 1
        for BulkIndex = 1: BulkFaultcount
            if Input_Bulk[BulkIndex,8] == -1.0
                Input_Bulk[BulkIndex,8] = 0.0
            else        
                Input_Bulk[BulkIndex,8] = 180.0
            end    
        end
    elseif Switch_StrikeSlip_or_ReverseNormal ==2
        for BulkIndex = 1: BulkFaultcount
            if Input_Bulk[BulkIndex,7] < 90.0
                if Input_Bulk[BulkIndex,8] == -1.0
                    Input_Bulk[BulkIndex,8] = 90.0
                else        
                    Input_Bulk[BulkIndex,8] = 270.0
                end    
            else 
                if Input_Bulk[BulkIndex,8] == -1.0
                    Input_Bulk[BulkIndex,8] = 270.0
                else        
                    Input_Bulk[BulkIndex,8] = 90.0
                end    
            end
        end
    end   
    
    return Input_Bulk

end


function ReadBulkInput(InputBulkFileName)

    ######################### Check if Rectangle or Triangle #############################
    RorT = ""
    Input_Bulk=readdlm(InputBulkFileName)
    if size(Input_Bulk, 2) == 18
        println("Rectangle")
        RorT = "R"
    elseif  size(Input_Bulk, 2) == 20
        println("Triangle")
        RorT = "T"
    else 
        error("Input Bulk Fault Geometry file should have 18 or 20 columns")
    end


    ########################## Read and Remove Header then Segmentize #########################
    Switch_StrikeSlip_or_ReverseNormal = Input_Bulk[2,1] 
    ShearModulus = Input_Bulk[2,2]
    PoissonRatio = Input_Bulk[2,3]
    RockDensity = Input_Bulk[2,4]
    DropCrit= Input_Bulk[2,5]
    DropCritNormalStressMultiplier= Input_Bulk[2,6]
    MinimumNS=Input_Bulk[2,7]
    Input_Bulk=Input_Bulk[4:end,:]
    LoadingFaultCount   = 0

    if RorT == "R"
        Input_Bulk=Input_Bulk[sortperm(Input_Bulk[:, 17]), :]  
        Input_Bulk = LRtoRake(Switch_StrikeSlip_or_ReverseNormal, Input_Bulk)# Adjust LRRN to rake angle    
        for i in eachindex(Input_Bulk[:,1]) # Warning if positive depth exists
            if Input_Bulk[i,3]<Input_Bulk[i,5]/2*sind(Input_Bulk[i,7])
                println("Caution! Fault ",i," may have negative depth")
            end
        end
        Input_Segment = BulkToSegment(Input_Bulk);
        FaultCount=   size(Input_Segment,1)
        Input_Segment = Input_Segment[sortperm(Input_Segment[:, 16]), :] # move the loading faults to the top
        LoadingFaultCount = sum(Input_Segment[:,16] .> 0)

    elseif RorT == "T"
        Input_Bulk = Input_Bulk[sortperm(Input_Bulk[:, 19]), :]     
        Input_Segment = Input_Bulk
        FaultCount=   size(Input_Segment,1)
        Input_Segment = Input_Segment[sortperm(Input_Segment[:, 19]), :] # move the loading faults to the top
        LoadingFaultCount = sum(Input_Segment[:,19] .> 0)

    end

    return Input_Segment, LoadingFaultCount, ShearModulus, PoissonRatio, RockDensity, 
           Switch_StrikeSlip_or_ReverseNormal, DropCrit, DropCritNormalStressMultiplier, MinimumNS, RorT, FaultCount

end



function ReadSegmentInput(Input_Segment, FaultCount, RorT)

    if RorT == "R"

        FaultCenter = Input_Segment[:,1:3]
        FaultLengthStrike = Input_Segment[:,4]
        FaultLengthDip = Input_Segment[:,5]
        FaultStrikeAngle = Input_Segment[:,6]
        FaultDipAngle = Input_Segment[:,7]
        FaultRakeAngle = Input_Segment[:,8]
        Fault_a = Input_Segment[:,9]
        Fault_b = Input_Segment[:,10]
        Fault_Dc = Input_Segment[:,11]
        Fault_Theta_i = Input_Segment[:,12]
        Fault_V_i = Input_Segment[:,13]
        Fault_Friction_i = Input_Segment[:,14]
        Fault_NormalStress = Input_Segment[:,15]
        Fault_V_Const = Input_Segment[:,16]
        Fault_BulkIndex = Input_Segment[:,17]
        FaultLengthStrike_Bulk = Input_Segment[:,18]
        FaultLengthDip_Bulk = Input_Segment[:,19]
        NormalStiffnessZero = 0

    else
        P1 = Input_Segment[:,1:3]
        P2 = Input_Segment[:,4:6]
        P3 = Input_Segment[:,7:9]
        FaultRakeAngle = Input_Segment[:,10]
        FaultCenter = zeros(FaultCount,3)
        FaultCenter[:,1] = (P1[:,1] + P2[:,1] + P3[:,1]) /3
        FaultCenter[:,2] = (P1[:,2] + P2[:,2] + P3[:,2]) /3
        FaultCenter[:,3] = (P1[:,3] + P2[:,3] + P3[:,3]) /3
        
        Fault_a = Input_Segment[:,11]
        Fault_b = Input_Segment[:,12]
        Fault_Dc = Input_Segment[:,13]
        Fault_Theta_i = Input_Segment[:,14]
        Fault_V_i = Input_Segment[:,15]
        Fault_Friction_i = Input_Segment[:,16]
        Fault_NormalStress = Input_Segment[:,17] - Input_Segment[:,18] .*  FaultCenter[:,3]
        Fault_V_Const = Input_Segment[:,19]
        Fault_BulkIndex = FaultRakeAngle * 0
        FaultLengthStrike = maximum.(eachrow(abs.([P1-P2 P2-P3 P1-P3])))
        FaultLengthDip = FaultLengthStrike
        FaultStrikeAngle = zeros(FaultCount)
        FaultDipAngle = zeros(FaultCount)
        FaultLengthStrike_Bulk = FaultLengthStrike
        FaultLengthDip_Bulk = FaultLengthDip
        NormalStiffnessZero = 0
    end

    return FaultCenter, Fault_a, Fault_b, Fault_Dc, Fault_Theta_i, Fault_V_i, Fault_Friction_i, Fault_NormalStress, 
            Fault_V_Const, Fault_BulkIndex, FaultLengthStrike, FaultLengthDip, FaultStrikeAngle, 
            FaultDipAngle, FaultRakeAngle, FaultLengthStrike_Bulk, FaultLengthDip_Bulk, NormalStiffnessZero
end





function RotVerts_UnitVectors(Input_Segment, FaultCount, Rake)

    P1 = Input_Segment[:,1:3]
    P2 = Input_Segment[:,4:6]
    P3 = Input_Segment[:,7:9]
    UnitVector_Normal = zeros(FaultCount,3)
    UnitVector_StrikeSlip = zeros(FaultCount,3)
    UnitVector_DipSlip = zeros(FaultCount,3)
    UnitVector_Slip = zeros(FaultCount,3)

        VertsReversed = 0
    for ElemIdx = 1:FaultCount
        P1_i = P1[ElemIdx,:] 
        P2_i = P2[ElemIdx,:] 
        P3_i = P3[ElemIdx,:] 

        UnitVector_Normal_i = cross(P2_i-P1_i, P3_i-P1_i) / norm(cross(P2_i-P1_i, P3_i-P1_i))

        if UnitVector_Normal_i[2] <= 0
            VertsReversed = VertsReversed +1 
            P_temp = P1_i
            P1_i = P2_i
            P2_i = P_temp        
        end
        P1[ElemIdx,:] = P1_i
        P2[ElemIdx,:] = P2_i
        P3[ElemIdx,:] = P3_i
        
        UnitVector_Normal[ElemIdx,:]= cross(P2_i-P1_i, P3_i-P1_i) / norm(cross(P2_i-P1_i, P3_i-P1_i))
        UnitVector_StrikeSlip[ElemIdx,:] = cross(UnitVector_Normal[ElemIdx,:], [0, 0, 1]) / 
                                            norm(cross(UnitVector_Normal[ElemIdx,:], [0, 0, 1]) )
        UnitVector_DipSlip[ElemIdx,:] = cross(UnitVector_StrikeSlip[ElemIdx,:], UnitVector_Normal[ElemIdx,:]) /
                                            norm(cross(UnitVector_StrikeSlip[ElemIdx,:], UnitVector_Normal[ElemIdx,:]))
        UnitVector_Slip[ElemIdx,:] = UnitVector_StrikeSlip[ElemIdx,:] * cosd(Rake[ElemIdx]) + 
                                        UnitVector_DipSlip[ElemIdx,:] * sind(Rake[ElemIdx])

    end
            println("Verts order reversed: ", VertsReversed)
    return P1, P2, P3, UnitVector_Normal, UnitVector_StrikeSlip, UnitVector_DipSlip, UnitVector_Slip
end




###################################################################
####### Build Stiffness Matrix  Stirke Slip Vector By Part  #######
###################################################################

function StiffnessMatrix_ByParts_Calculation_Rec(Input_SegmentSource, Input_SegmentReceiver, ShearModulus, PoissonRatio,
    CurrentPart, TotalParts)
    FaultCountSource=size(Input_SegmentSource,1)
    FaultCenterSource=Input_SegmentSource[:,1:3]
    FaultLengthStrikeSource=Input_SegmentSource[:,4]
    FaultLengthDipSource=Input_SegmentSource[:,5]
    FaultStrikeAngleSource=Input_SegmentSource[:,6]
    FaultDipAngleSource=Input_SegmentSource[:,7]
    FaultRakeSource=Input_SegmentSource[:,8]

    FaultCountReceiver=size(Input_SegmentReceiver,1)
    FaultCenterReceiver=Input_SegmentReceiver[:,1:3]
    FaultLengthStrikeReceiver=Input_SegmentReceiver[:,4]
    FaultLengthDipReceiver=Input_SegmentReceiver[:,5]
    FaultStrikeAngleReceiver=Input_SegmentReceiver[:,6]
    FaultDipAngleReceiver=Input_SegmentReceiver[:,7]
    FaultRakeReceiver=Input_SegmentReceiver[:,8]
    # println(FaultCountSource, "  ", FaultCountReceiver)

    StiffnessMatrixShear = zeros(FaultCountReceiver,FaultCountSource)
    StiffnessMatrixNormal = zeros(FaultCountReceiver,FaultCountSource)
    # println(CurrentPart,"/",TotalParts)
    
        for SourceIndex=1:FaultCountSource;
            
            # println(SourceIndex,"  ",CurrentPart,"/",TotalParts)

            ####################################
            ##### get source geometry and slip

            SourceCenter = FaultCenterSource[SourceIndex,:];
            SourceLengthStrike = FaultLengthStrikeSource[SourceIndex];
            SourceLengthDip = FaultLengthDipSource[SourceIndex];
            SourceStrikeAngle = FaultStrikeAngleSource[SourceIndex];
            SourceDipAngle = FaultDipAngleSource[SourceIndex];
            SourceRake = FaultRakeSource[SourceIndex];
                    
            ReceiverCenter = FaultCenterReceiver;
            ReceiverStrikeAngle = FaultStrikeAngleReceiver;
            ReceiverDipAngle = FaultDipAngleReceiver;
            ReceiverRake = FaultRakeReceiver;
            RelativeStrkieAngle = ReceiverStrikeAngle .- SourceStrikeAngle
            

            DISL1 = cosd(SourceRake); # Left Latteral is +1 for Okada
            DISL2 = sind(SourceRake);
            DISL3 = 0;            
                    
            DEPTH=SourceCenter[3]; # Source Depth
            AL1=SourceLengthStrike/2;
            AL2=SourceLengthStrike/2;
            AW1=SourceLengthDip/2;
            AW2=SourceLengthDip/2;
            
            LameFirstParam=2*ShearModulus*PoissonRatio/(1-2*PoissonRatio);
            ALPHA=(LameFirstParam+ShearModulus)/(LameFirstParam+2*ShearModulus);
            

            #######################################################################
            ##### Calculate Receiver Point Relative to the Source and Source frame
            
            X_Dist = ReceiverCenter[:,1] .- SourceCenter[1];
            Y_Dist = ReceiverCenter[:,2] .- SourceCenter[2];
            
            X = X_Dist .* cosd(-SourceStrikeAngle) .- Y_Dist .* sind(-SourceStrikeAngle);
            Y = X_Dist .* sind(-SourceStrikeAngle) .+ Y_Dist .* cosd(-SourceStrikeAngle);
            Z = -ReceiverCenter[:,3];


            #######################################################################
            ##### Calculate Stress Change at Source Frame

            UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET = Okada_DC3D_Vector(ALPHA,
                X,Y,Z,DEPTH,SourceDipAngle,
                AL1,AL2,AW1,AW2,DISL1,DISL2,DISL3);

            StressXX_SourceFrame=(LameFirstParam*(UXX+UYY+UZZ) + 2*ShearModulus*UXX);
            StressYY_SourceFrame=(LameFirstParam*(UXX+UYY+UZZ) + 2*ShearModulus*UYY);
            StressZZ_SourceFrame=(LameFirstParam*(UXX+UYY+UZZ) + 2*ShearModulus*UZZ);
            StressXY_SourceFrame=(UXY + UYX)*ShearModulus;
            StressXZ_SourceFrame=(UXZ + UZX)*ShearModulus;
            StressYZ_SourceFrame=(UYZ + UZY)*ShearModulus;     
            


            for ReceiverIdx = 1:FaultCountReceiver
                Stress_SourceFrame=[StressXX_SourceFrame[ReceiverIdx] StressXY_SourceFrame[ReceiverIdx] StressXZ_SourceFrame[ReceiverIdx]
                StressXY_SourceFrame[ReceiverIdx] StressYY_SourceFrame[ReceiverIdx] StressYZ_SourceFrame[ReceiverIdx]
                StressXZ_SourceFrame[ReceiverIdx] StressYZ_SourceFrame[ReceiverIdx] StressZZ_SourceFrame[ReceiverIdx]];
        

                #######################################################################
                ##### Rotate Source Frame Stress to Flat Receiver 

                RotationMat_FromReceiver_Strike=
                [cosd(-RelativeStrkieAngle[ReceiverIdx]) -sind(-RelativeStrkieAngle[ReceiverIdx])  0
                sind(-RelativeStrkieAngle[ReceiverIdx]) cosd(-RelativeStrkieAngle[ReceiverIdx]) 0
                0  0  1];

                RotationMat_FromReceiver_Dip=
                [1 0 0
                0 cosd(-ReceiverDipAngle[ReceiverIdx]) -sind(-ReceiverDipAngle[ReceiverIdx])
                0 sind(-ReceiverDipAngle[ReceiverIdx]) cosd(-ReceiverDipAngle[ReceiverIdx])]
                
                RotationMat_FromReceiver_All = RotationMat_FromReceiver_Dip*RotationMat_FromReceiver_Strike
                                
                Stress_Receiver = RotationMat_FromReceiver_All*Stress_SourceFrame*RotationMat_FromReceiver_All'

                #######################################################################
                ##### Read Normal and Shear Stress Change
                StiffnessMatrixNormal[ReceiverIdx,SourceIndex] = - Stress_Receiver[3,3]  # compression is negative
                StiffnessMatrixShear[ReceiverIdx,SourceIndex] = cosd(ReceiverRake[ReceiverIdx]) * Stress_Receiver[1,3] + sind(ReceiverRake[ReceiverIdx]) * Stress_Receiver[2,3]  # right latteral is negative

                # println(SourceDipAngle,"  ",Z,"  ",DEPTH, " ", StressZZ_SourceFrame, "  ",Stress_Receiver[3,3])
            end
        
        if FaultCountSource > 500
            if SourceIndex % 50 == 0
                println("Completed Source Element: ", SourceIndex, "/", FaultCountSource)
            end
        end
        end
        print("\033c")                  


    return StiffnessMatrixShear, StiffnessMatrixNormal 
end
 


function BuildMatrixByParts_Even_Rec(FaultCount, Input_Segment,  ShearModulus, PoissonRatio)
    ElementPartRoughCount = 2000
    DivisionCount = round(Int,FaultCount / ElementPartRoughCount)
    if DivisionCount == 0; DivisionCount =1; end
    PartedElementCount = FaultCount ÷ DivisionCount
    StiffnessMatrixShearOriginal= zeros(FaultCount,FaultCount)
    StiffnessMatrixNormalOriginal = zeros(FaultCount,FaultCount)
    TotalParts = DivisionCount^2
    CurrentPart = 0
    println("preparing for discretization by parts. Total Parts ", TotalParts)
    for i=1:DivisionCount
        for j=1:DivisionCount
            

            CurrentPart =  CurrentPart +1
            
            println("Part: ", CurrentPart,"/",DivisionCount*DivisionCount)

            Init_S = (i-1)*PartedElementCount+1
            Fin_S = i*PartedElementCount
            Init_R =  (j-1)*PartedElementCount+1
            Fin_R = j*PartedElementCount
            if i == DivisionCount; Fin_S = FaultCount; end
            if j == DivisionCount; Fin_R = FaultCount; end
    
            StiffnessMatrixShearOriginal[Init_R:Fin_R,Init_S:Fin_S], StiffnessMatrixNormalOriginal[Init_R:Fin_R,Init_S:Fin_S] = 
            StiffnessMatrix_ByParts_Calculation_Rec(Input_Segment[Init_S:Fin_S,:], Input_Segment[Init_R:Fin_R,:], ShearModulus, PoissonRatio,
                                                CurrentPart, TotalParts)
        end
    end
    
    return StiffnessMatrixShearOriginal, StiffnessMatrixNormalOriginal
end


function StiffnessMatrix_ByParts_Calculation_Tri(P1, P2, P3, Rake, FaultCenter, ShearModulus, lambda,
                            UnitVector_Normal, UnitVector_Slip)
    

    FaultCountSource=size(P1,1)
    FaultCountReceiver=size(FaultCenter,1)
    StiffnessMatrix_Normal_Part= zeros(FaultCountReceiver, FaultCountSource)
    StiffnessMatrix_Shear_Part = zeros(FaultCountReceiver, FaultCountSource)
    for ElemIdx = 1:FaultCountSource
    p1 = P1[ElemIdx,:]
    p2 = P2[ElemIdx,:]
    p3 = P3[ElemIdx,:]
    Ss = -cosd(Rake[ElemIdx])
    Ds = sind(Rake[ElemIdx])
    # println(Ss, "  ", Ds)
    Ts = 0.0
        Stress,Strain = TDstressHS(FaultCenter[:,1],FaultCenter[:,2],FaultCenter[:,3],p1,p2,p3,
        Ss,Ds,Ts,ShearModulus,lambda)
        # println(ElemIdx)   
        for ElemIdx2 = 1:FaultCountReceiver
            Stress_i = [Stress[ElemIdx2,1] Stress[ElemIdx2,4] Stress[ElemIdx2,5]
                        Stress[ElemIdx2,4] Stress[ElemIdx2,2] Stress[ElemIdx2,6]
                        Stress[ElemIdx2,5] Stress[ElemIdx2,6] Stress[ElemIdx2,3]]

            TVector = Stress_i * UnitVector_Normal[ElemIdx2,:]
            Stress_Normal = dot(TVector, UnitVector_Normal[ElemIdx2,:])
            # Stress_SS = dot(TVector, UnitVector_StrikeSlip[ElemIdx2,:])
            # Stress_Dip = dot(TVector, UnitVector_DipSlip[ElemIdx2,:])
            Stress_Shear = dot(TVector, UnitVector_Slip[ElemIdx2,:])
            StiffnessMatrix_Normal_Part[ElemIdx2, ElemIdx] = Stress_Normal
            StiffnessMatrix_Shear_Part[ElemIdx2, ElemIdx] = Stress_Shear
        end

        if FaultCountSource > 500
            if ElemIdx % 50 == 0
                println("Completed Source Element: ", ElemIdx, "/", FaultCountSource)
            end
        end
    end
    # println("Fault Count Source: ", FaultCountSource, " Fault Count Receiver: ", FaultCountReceiver)
        print("\033c")                  
    return StiffnessMatrix_Shear_Part, StiffnessMatrix_Normal_Part

end 




function BuildMatrixByParts_Even_Tri(P1, P2, P3, Rake, FaultCenter, UnitVector_Normal, UnitVector_Slip,
                                FaultCount, ShearModulus, PoissonRatio)

    ElementPartRoughCount = 2000

    StiffnessMatrix_Shear=zeros(FaultCount,FaultCount)
    StiffnessMatrix_Normal=zeros(FaultCount,FaultCount)

    DivisionCount = round(Int,FaultCount / ElementPartRoughCount)
    if DivisionCount == 0; DivisionCount =1; end
    PartedElementCount = FaultCount ÷ DivisionCount
    TotalParts = DivisionCount^2
    CurrentPart = 0
    lambda = 2 * ShearModulus * PoissonRatio / (1 - 2 * PoissonRatio)
    println("preparing for discretization by parts. Total Parts ", TotalParts)
    println("Compiling Stiffness Matrix Function. This may take a while if first run")

    for i=1:DivisionCount
        for j=1:DivisionCount
            CurrentPart =  CurrentPart +1 
            
            println("Part: ", CurrentPart,"/",DivisionCount*DivisionCount)
            Init_S = (i-1)*PartedElementCount+1
            Fin_S = i*PartedElementCount
            Init_R =  (j-1)*PartedElementCount+1
            Fin_R = j*PartedElementCount
            if i == DivisionCount; Fin_S = FaultCount; end
            if j == DivisionCount; Fin_R = FaultCount; end

            StiffnessMatrix_Shear[Init_R:Fin_R,Init_S:Fin_S], StiffnessMatrix_Normal[Init_R:Fin_R,Init_S:Fin_S] = 
            StiffnessMatrix_ByParts_Calculation_Tri(P1[Init_S:Fin_S,:], P2[Init_S:Fin_S,:], P3[Init_S:Fin_S,:], Rake[Init_S:Fin_S,:],
                                                    FaultCenter[Init_R:Fin_R,:], ShearModulus, lambda,  
                                                    UnitVector_Normal[Init_R:Fin_R,:], UnitVector_Slip[Init_R:Fin_R,:])    
                                        
        end
    end

    
    return StiffnessMatrix_Shear, StiffnessMatrix_Normal
end




function HmatBuild_T(ShearModulus, PoissonRatio, ElementRange_SR, FaultRakeAngle, FaultCenter,
                        UnitVector_Normal, UnitVector_Slip, Admissible, Tolerance, P1, P2, P3)

    ################################  Discritize ##############################
    ElementPartRoughCount = 2000

    lambda = 2 * ShearModulus * PoissonRatio / (1 - 2 * PoissonRatio)
    BlockCount = length(ElementRange_SR[:,1])
    ShearStiffness_H = Any[0]
    NormalStiffness_H = Any[0]
    Ranks_Shear = zeros(Int, BlockCount)
    Ranks_Normal = zeros(Int, BlockCount)
    FaultCount = size(FaultCenter,1)
    TotalElments = FaultCount * FaultCount 
    println("Building Hmatrix Block by Block")
    println("Full Matrix Will not be saved")
    println("Preparing for discretization")
    println("Initiation of first part discretization may take a few minutes due to compliation")
    BlockSize = 0.0

    for BlockIndex = 1: BlockCount
        P1_S = P1[ElementRange_SR[BlockIndex,1]:ElementRange_SR[BlockIndex,2],:]
        P2_S = P2[ElementRange_SR[BlockIndex,1]:ElementRange_SR[BlockIndex,2],:]
        P3_S = P3[ElementRange_SR[BlockIndex,1]:ElementRange_SR[BlockIndex,2],:]
        FaultRakeAngle_S = FaultRakeAngle[ElementRange_SR[BlockIndex,1]:ElementRange_SR[BlockIndex,2]]
        FaultCenter_R = FaultCenter[ElementRange_SR[BlockIndex,3]:ElementRange_SR[BlockIndex,4],:]
        UnitVector_Normal_R = UnitVector_Normal[ElementRange_SR[BlockIndex,3]:ElementRange_SR[BlockIndex,4],:]
        UnitVector_Slip_R = UnitVector_Slip[ElementRange_SR[BlockIndex,3]:ElementRange_SR[BlockIndex,4],:]

        SourceCount = ElementRange_SR[BlockIndex,2] - ElementRange_SR[BlockIndex,1] + 1
        ReceiverCount = ElementRange_SR[BlockIndex,4] - ElementRange_SR[BlockIndex,3] + 1
        StiffnessMatrixShearThisBlock=zeros(ReceiverCount,SourceCount)
        StiffnessMatrixNormalThisBlock=zeros(ReceiverCount,SourceCount)
        DivisionCountS = round(Int,SourceCount / ElementPartRoughCount)
        DivisionCountR = round(Int,ReceiverCount / ElementPartRoughCount)
        if DivisionCountS == 0; DivisionCountS =1; end
        if DivisionCountR == 0; DivisionCountR =1; end
        PartedElementCountS = SourceCount ÷ DivisionCountS
        PartedElementCountR = ReceiverCount ÷ DivisionCountR
        TotalParts = DivisionCountS * DivisionCountR
        CurrentPart = 0
        for i=1:DivisionCountS
            for j=1:DivisionCountR
                CurrentPart =  CurrentPart +1
                Init_S = (i-1)*PartedElementCountS + 1
                Fin_S = i*PartedElementCountS
                Init_R =  (j-1)*PartedElementCountR + 1
                Fin_R = j*PartedElementCountR
                if i == DivisionCountS; Fin_S = SourceCount; end
                if j == DivisionCountR; Fin_R = ReceiverCount; end
                BlockSize = BlockSize + (Fin_S - Init_S) * (Fin_R - Init_R)
                println("Part: ", CurrentPart,"/",DivisionCountS*DivisionCountR, " BlockIndex: ",BlockIndex, "/",BlockCount," Progress:",BlockSize/TotalElments)

                StiffnessMatrixShearThisBlock[Init_R:Fin_R,Init_S:Fin_S], StiffnessMatrixNormalThisBlock[Init_R:Fin_R,Init_S:Fin_S] = 
                StiffnessMatrix_ByParts_Calculation_Tri(P1_S[Init_S:Fin_S,:], P2_S[Init_S:Fin_S,:], P3_S[Init_S:Fin_S,:], FaultRakeAngle_S[Init_S:Fin_S,:],
                                                        FaultCenter_R[Init_R:Fin_R,:], ShearModulus, lambda,  
                                                        UnitVector_Normal_R[Init_R:Fin_R,:], UnitVector_Slip_R[Init_R:Fin_R,:])   
            end
        end                        

            ######################################################################
            ############################# Compress ###############################
        if Admissible[BlockIndex] > 0
            ApproxMatrixS = pqrfact(StiffnessMatrixShearThisBlock, atol = Tolerance)
            push!(ShearStiffness_H,ApproxMatrixS)
            Ranks_Shear[BlockIndex] = size(ApproxMatrixS[:Q],2)
            
            ApproxMatrixN = pqrfact(StiffnessMatrixNormalThisBlock, atol = Tolerance)
            push!(NormalStiffness_H,ApproxMatrixN)
            Ranks_Normal[BlockIndex] = size(ApproxMatrixN[:Q],2)
        else 
            push!(ShearStiffness_H,StiffnessMatrixShearThisBlock)
            push!(NormalStiffness_H,StiffnessMatrixNormalThisBlock)
        end
    end
    ShearStiffness_H = ShearStiffness_H[2:end]
    NormalStiffness_H = NormalStiffness_H[2:end]

    return ShearStiffness_H, NormalStiffness_H, Ranks_Shear, Ranks_Normal
end


function HmatBuild_R(ShearModulus, PoissonRatio, ElementRange_SR,Input_Segment, Admissible, Tolerance)

    FaultCount = size(Input_Segment)[1]
    ################################  Discritize ##############################
    ElementPartRoughCount = 2000
    BlockCount = length(ElementRange_SR[:,1])
    ShearStiffness_H = Any[0]
    NormalStiffness_H = Any[0]
    Ranks_Shear = zeros(Int, BlockCount)
    Ranks_Normal = zeros(Int, BlockCount)

    TotalElments = FaultCount * FaultCount 
    println("Building Hmatrix Block by Block")
    println("Full Matrix Will not be saved")
    println("Preparing for discretization")
    BlockSize = 0.0



    for BlockIndex = 1: BlockCount

        Input_SegmentS = Input_Segment[ElementRange_SR[BlockIndex,1]:ElementRange_SR[BlockIndex,2],:]
        Input_SegmentR = Input_Segment[ElementRange_SR[BlockIndex,3]:ElementRange_SR[BlockIndex,4],:]

        SourceCount = ElementRange_SR[BlockIndex,2] - ElementRange_SR[BlockIndex,1] + 1
        ReceiverCount = ElementRange_SR[BlockIndex,4] - ElementRange_SR[BlockIndex,3] + 1
        StiffnessMatrixShearThisBlock=zeros(ReceiverCount,SourceCount)
        StiffnessMatrixNormalThisBlock=zeros(ReceiverCount,SourceCount)
        DivisionCountS = round(Int,SourceCount / ElementPartRoughCount)
        DivisionCountR = round(Int,ReceiverCount / ElementPartRoughCount)
        if DivisionCountS == 0; DivisionCountS =1; end
        if DivisionCountR == 0; DivisionCountR =1; end
        PartedElementCountS = SourceCount ÷ DivisionCountS
        PartedElementCountR = ReceiverCount ÷ DivisionCountR
        TotalParts = DivisionCountS * DivisionCountR
        CurrentPart = 0
        for i=1:DivisionCountS
            for j=1:DivisionCountR
                CurrentPart =  CurrentPart +1
                Init_S = (i-1)*PartedElementCountS + 1
                Fin_S = i*PartedElementCountS
                Init_R =  (j-1)*PartedElementCountR + 1
                Fin_R = j*PartedElementCountR
                if i == DivisionCountS; Fin_S = SourceCount; end
                if j == DivisionCountR; Fin_R = ReceiverCount; end
                BlockSize = BlockSize + (Fin_S - Init_S) * (Fin_R - Init_R)
                println("Part: ", CurrentPart,"/",DivisionCountS*DivisionCountR, " BlockIndex: ",BlockIndex, "/",BlockCount," Progress:",BlockSize/TotalElments)

                StiffnessMatrixShearThisBlock[Init_R:Fin_R,Init_S:Fin_S], StiffnessMatrixNormalThisBlock[Init_R:Fin_R,Init_S:Fin_S] = 
                StiffnessMatrix_ByParts_Calculation_Rec(Input_SegmentS[Init_S:Fin_S,:], Input_SegmentR[Init_R:Fin_R,:], ShearModulus, PoissonRatio,
                                                    CurrentPart, TotalParts)

                                          
            
            end
        end                        

            ######################################################################
            ############################# Compress ###############################
        if Admissible[BlockIndex] > 0
            ApproxMatrixS = pqrfact(StiffnessMatrixShearThisBlock, atol = Tolerance)
            push!(ShearStiffness_H,ApproxMatrixS)
            Ranks_Shear[BlockIndex] = size(ApproxMatrixS[:Q],2)
            
            ApproxMatrixN = pqrfact(StiffnessMatrixNormalThisBlock, atol = Tolerance)
            push!(NormalStiffness_H,ApproxMatrixN)
            Ranks_Normal[BlockIndex] = size(ApproxMatrixN[:Q],2)
        else 
            push!(ShearStiffness_H,StiffnessMatrixShearThisBlock)
            push!(NormalStiffness_H,StiffnessMatrixNormalThisBlock)
        end
    # ApproxMatrixS*ones(200)
    end
                ShearStiffness_H = ShearStiffness_H[2:end]
                NormalStiffness_H = NormalStiffness_H[2:end]


    return  ShearStiffness_H, NormalStiffness_H, Ranks_Shear, Ranks_Normal
end