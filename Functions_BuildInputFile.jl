

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


###################################################################
############### Check Orientation with Loading Faults #############
################## This is Currently Not being Used ###############
###################################################################
   
function Function_CheckOrientation(Input_Segment, ShearModulus, PoissonRatio)
    ## InputFile Order 
    ##  1.Ctr_X     2.Ctr_Y 3.Ctr_Z 4.St_L	    5.Dip_L	    6.StAng	    7.DipAng	8.LR/RN
    ##  9.a         10.b	11.Dc	12.Theta_i	13. V_i     14. Friction_i 15.NormalStress  
    ##  16. V_Const 17. Bulk Number     18. Bulk Strike Length      19. Bulk Dip Length
    
    FaultCount=size(Input_Segment,1)
    Fault_V_Const=Input_Segment[:,16]
    LoadingFaultCount=length(Fault_V_Const[Fault_V_Const.>0])
    ReceiverOrientation=zeros(FaultCount,LoadingFaultCount)
    StiffnessMatrixNormal=zeros(FaultCount,FaultCount)
    StiffnessMatrixShear=zeros(FaultCount,FaultCount)

    FaultCenter=Input_Segment[:,1:3]
    FaultLengthStrike=Input_Segment[:,4]
    FaultLengthDip=Input_Segment[:,5]
    FaultStrikeAngle=Input_Segment[:,6] 
    FaultDipAngle=Input_Segment[:,7]
    FaultLLRR=Input_Segment[:,8]

    println("Checking Fault Orientations")

    for ReceiverIndex=1:FaultCount;
    
        SourceNumber=0
        for SourceIndex=FaultCount-LoadingFaultCount+1:FaultCount;
            SourceNumber=SourceNumber+1;
            
            
            SourceCenter=FaultCenter[SourceIndex,:];
            SourceLengthStrike=FaultLengthStrike[SourceIndex];
            SourceLengthDip=FaultLengthDip[SourceIndex];
            SourceStrikeAngle=FaultStrikeAngle[SourceIndex];
            SourceDipAngle=FaultDipAngle[SourceIndex];
            SourceLLRR=FaultLLRR[SourceIndex];
                    
            ReceiverCenter=FaultCenter[ReceiverIndex,:];
            ReceiverLengthStrike=FaultLengthStrike[ReceiverIndex];
            ReceiverLengthDip=FaultLengthDip[ReceiverIndex];
            ReceiverStrikeAngle=FaultStrikeAngle[ReceiverIndex];
            ReceiverDipAngle=FaultDipAngle[ReceiverIndex];
            ReceiverLLRR=FaultLLRR[ReceiverIndex];
            DISL1=-SourceLLRR;
            DISL2=0;
            DISL3=0;
                    
            Z=-ReceiverCenter[3]; # Observation Depth
            DEPTH=SourceCenter[3];
            AL1=SourceLengthStrike/2;
            AL2=SourceLengthStrike/2;
            AW1=SourceLengthDip/2;
            AW2=SourceLengthDip/2;
            
            LameFirstParam=2*ShearModulus*PoissonRatio/(1-2*PoissonRatio);
            ALPHA=(LameFirstParam+ShearModulus)/(LameFirstParam+2*ShearModulus);
            
            # Rotation Matrix Source
            RotMat_Source_Strike=[cosd(SourceStrikeAngle) -sind(SourceStrikeAngle) 0
                sind(SourceStrikeAngle) cosd(SourceStrikeAngle) 0
                0  0  1];
            RotMat_Source_Dip=[1  0  0
                0  cosd(SourceDipAngle) -sind(SourceDipAngle)
                0  sind(SourceDipAngle) cosd(SourceDipAngle)];
            
            # Rotation Matrix Receiver
            RotationMat_FromReceiver_Strike=
            [cosd(-ReceiverStrikeAngle) -sind(-ReceiverStrikeAngle)  0
            sind(-ReceiverStrikeAngle) cosd(-ReceiverStrikeAngle) 0
            0  0  1];
            RotationMat_FromReceiver_Dip=
            [1 0 0
            0 cosd(-ReceiverDipAngle) -sind(-ReceiverDipAngle)
            0 sind(-ReceiverDipAngle) cosd(-ReceiverDipAngle)]
            RotationMat_FromReceiver_All=RotationMat_FromReceiver_Dip*RotationMat_FromReceiver_Strike;
            
            
            # Point Calculation
            
            X_Receiver=ReceiverCenter[1];
            Y_Receiver=ReceiverCenter[2];
            X_Dist=X_Receiver-SourceCenter[1];
            Y_Dist=Y_Receiver-SourceCenter[2];
            
            X=X_Dist*cosd(-SourceStrikeAngle)-Y_Dist*sind(-SourceStrikeAngle);
            Y=X_Dist*sind(-SourceStrikeAngle)+Y_Dist*cosd(-SourceStrikeAngle);
            Z=-ReceiverCenter[3];
            
            UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET = Okada_DC3D(ALPHA,
                X,Y,Z,DEPTH,SourceDipAngle,
                AL1,AL2,AW1,AW2,DISL1,DISL2,DISL3);
            
            UnRotatedGradMat=
            [UXX  UYX  UZX
            UXY  UYY  UZY
            UXZ  UYZ  UZZ];
            #DispRefAxisRotatd=RotMat_Source_Strike*[UX;UY;UZ];
            GradientRefAxisOrigin=RotMat_Source_Strike*UnRotatedGradMat*RotMat_Source_Strike';
            
            #Result_UX=DispRefAxisRotatd[1];
            #Result_UY=DispRefAxisRotatd[2];
            #Result_UZ=DispRefAxisRotatd[3];
            Result_UXX=GradientRefAxisOrigin[1,1];
            Result_UYX=GradientRefAxisOrigin[2,1];
            Result_UZX=GradientRefAxisOrigin[3,1];
            Result_UXY=GradientRefAxisOrigin[1,2];
            Result_UYY=GradientRefAxisOrigin[2,2];
            Result_UZY=GradientRefAxisOrigin[3,2];
            Result_UXZ=GradientRefAxisOrigin[1,3];
            Result_UYZ=GradientRefAxisOrigin[2,3];
            Result_UZZ=GradientRefAxisOrigin[3,3];
            #Result_XLoc=X_Receiver;
            #Result_YLoc=Y_Receiver;
            #Result_ZLoc=Z;
            
            
            StressXX_SourceFrame=(LameFirstParam*(UXX+UYY+UZZ) + 2*ShearModulus*UXX)
            StressYY_SourceFrame=(LameFirstParam*(UXX+UYY+UZZ) + 2*ShearModulus*UYY)
            StressZZ_SourceFrame=(LameFirstParam*(UXX+UYY+UZZ) + 2*ShearModulus*UZZ)
            StressXY_SourceFrame=(Result_UXY + Result_UYX)*ShearModulus;
            StressXZ_SourceFrame=(Result_UXZ + Result_UZX)*ShearModulus;
            StressYZ_SourceFrame=(Result_UYZ + Result_UZY)*ShearModulus;
            
            Stress_SourceFrame=[StressXX_SourceFrame StressXY_SourceFrame StressXZ_SourceFrame
                StressXY_SourceFrame StressYY_SourceFrame StressYZ_SourceFrame
                StressXZ_SourceFrame StressYZ_SourceFrame StressZZ_SourceFrame];
            
            Stress_Receiver=RotationMat_FromReceiver_All*Stress_SourceFrame*RotationMat_FromReceiver_All';
            StiffnessMatrixNormal[ReceiverIndex,SourceIndex] = -Stress_Receiver[3,3];
            StiffnessMatrixShear[ReceiverIndex,SourceIndex] = ReceiverLLRR * Stress_Receiver[1,3];
            
            ReceiverOrientation[ReceiverIndex, SourceNumber]=StiffnessMatrixShear[ReceiverIndex,SourceIndex]/abs(StiffnessMatrixShear[ReceiverIndex,SourceIndex])
        end
        
    end

    return ReceiverOrientation

end



###################################################################
####### Build Stiffness Matrix  Stirke Slip Vector By Part  #######
###################################################################

function StiffnessMatrix_ByParts_Calculation_StrikeSlip(Input_SegmentSource, Input_SegmentReceiver, ShearModulus, PoissonRatio,
    CurrentPart, TotalParts)
    FaultCountSource=size(Input_SegmentSource,1)
    FaultCenterSource=Input_SegmentSource[:,1:3]
    FaultLengthStrikeSource=Input_SegmentSource[:,4]
    FaultLengthDipSource=Input_SegmentSource[:,5]
    FaultStrikeAngleSource=Input_SegmentSource[:,6]
    FaultDipAngleSource=Input_SegmentSource[:,7]
    FaultLLRRSource=Input_SegmentSource[:,8]

    FaultCountReceiver=size(Input_SegmentReceiver,1)
    FaultCenterReceiver=Input_SegmentReceiver[:,1:3]
    FaultLengthStrikeReceiver=Input_SegmentReceiver[:,4]
    FaultLengthDipReceiver=Input_SegmentReceiver[:,5]
    FaultStrikeAngleReceiver=Input_SegmentReceiver[:,6]
    FaultDipAngleReceiver=Input_SegmentReceiver[:,7]
    FaultLLRRReceiver=Input_SegmentReceiver[:,8]
    # println(FaultCountSource, "  ", FaultCountReceiver)

    StiffnessMatrixShear = zeros(FaultCountReceiver,FaultCountSource)
    StiffnessMatrixNormal = zeros(FaultCountReceiver,FaultCountSource)
    # println(CurrentPart,"/",TotalParts)
    
        for SourceIndex=1:FaultCountSource;
            
            println(SourceIndex,"  ",CurrentPart,"/",TotalParts)

            ####################################
            ##### get source geometry and slip

            SourceCenter = FaultCenterSource[SourceIndex,:];
            SourceLengthStrike = FaultLengthStrikeSource[SourceIndex];
            SourceLengthDip = FaultLengthDipSource[SourceIndex];
            SourceStrikeAngle = FaultStrikeAngleSource[SourceIndex];
            SourceDipAngle = FaultDipAngleSource[SourceIndex];
            SourceLLRR = FaultLLRRSource[SourceIndex];
                    
            ReceiverCenter = FaultCenterReceiver;
            ReceiverStrikeAngle = FaultStrikeAngleReceiver;
            ReceiverDipAngle = FaultDipAngleReceiver;
            ReceiverLLRR = FaultLLRRReceiver;
            RelativeStrkieAngle = ReceiverStrikeAngle .- SourceStrikeAngle
            

            DISL1 = -SourceLLRR; # Left Latteral is +1 for Okada
            DISL2 = 0;
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
                StiffnessMatrixShear[ReceiverIdx,SourceIndex] = - ReceiverLLRR[ReceiverIdx] * Stress_Receiver[1,3]  # right latteral is negative

                # println(SourceDipAngle,"  ",Z,"  ",DEPTH, " ", StressZZ_SourceFrame, "  ",Stress_Receiver[3,3])
            end
        
        end
        print("\033c")                  


    return StiffnessMatrixShear, StiffnessMatrixNormal 
end
 

###################################################################
####### Build Stiffness Matrix  Stirke Slip Vector By Part  #######
###################################################################

function StiffnessMatrix_ByParts_Calculation_NormalReverse(Input_SegmentSource, Input_SegmentReceiver, ShearModulus, PoissonRatio,
    CurrentPart, TotalParts)
    FaultCountSource=size(Input_SegmentSource,1)
    FaultCenterSource=Input_SegmentSource[:,1:3]
    FaultLengthStrikeSource=Input_SegmentSource[:,4]
    FaultLengthDipSource=Input_SegmentSource[:,5]
    FaultStrikeAngleSource=Input_SegmentSource[:,6]
    FaultDipAngleSource=Input_SegmentSource[:,7]
    FaultLLRRSource=Input_SegmentSource[:,8]

    FaultCountReceiver=size(Input_SegmentReceiver,1)
    FaultCenterReceiver=Input_SegmentReceiver[:,1:3]
    FaultLengthStrikeReceiver=Input_SegmentReceiver[:,4]
    FaultLengthDipReceiver=Input_SegmentReceiver[:,5]
    FaultStrikeAngleReceiver=Input_SegmentReceiver[:,6]
    FaultDipAngleReceiver=Input_SegmentReceiver[:,7]
    FaultLLRRReceiver=Input_SegmentReceiver[:,8]
    # println(FaultCountSource, "  ", FaultCountReceiver)


    StiffnessMatrixShear = zeros(FaultCountReceiver,FaultCountSource)
    StiffnessMatrixNormal = zeros(FaultCountReceiver,FaultCountSource)
    
    for SourceIndex=1:FaultCountSource;
        println(SourceIndex,"  ",CurrentPart,"/",TotalParts)
            

        ####################################
        ##### get source geometry and slip

        SourceCenter = FaultCenterSource[SourceIndex,:];
        SourceLengthStrike = FaultLengthStrikeSource[SourceIndex];
        SourceLengthDip = FaultLengthDipSource[SourceIndex];
        SourceStrikeAngle = FaultStrikeAngleSource[SourceIndex];
        SourceDipAngle = FaultDipAngleSource[SourceIndex];
        SourceLLRR = FaultLLRRSource[SourceIndex];
                
        ReceiverCenter = FaultCenterReceiver;
        ReceiverStrikeAngle = FaultStrikeAngleReceiver;
        ReceiverDipAngle = FaultDipAngleReceiver;
        ReceiverLLRR = FaultLLRRReceiver;
        RelativeStrkieAngle = ReceiverStrikeAngle .- SourceStrikeAngle
        
        
        if SourceDipAngle <= 90
            OrientationCoeff_Source = -1
        else
            OrientationCoeff_Source = 1
        end

        DISL1 = 0
        DISL2 = OrientationCoeff_Source * SourceLLRR # Reverse Slip is +1 for Okada
        DISL3 = 0                
                
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

            if ReceiverDipAngle[ReceiverIdx] <= 90
                OrientationCoeff_Receiver = -1
            else
                OrientationCoeff_Receiver = 1
            end

            StiffnessMatrixNormal[ReceiverIdx,SourceIndex] = - Stress_Receiver[3,3]  # compression is negative
            StiffnessMatrixShear[ReceiverIdx,SourceIndex] = OrientationCoeff_Receiver * ReceiverLLRR[ReceiverIdx] * Stress_Receiver[2,3]  # normal is negative

            # println(SourceDipAngle,"  ",Z,"  ",DEPTH, " ", StressZZ_SourceFrame, "  ",Stress_Receiver[3,3])
        end
    end                  
    # print("\033c")


return StiffnessMatrixShear, StiffnessMatrixNormal 
end
 



###################################################################
########### Remove Faults If too strongly interacting #############
###################################################################

function CheckTooClose(StiffnessMatrixShearOriginal, StiffnessMatrixNormalOriginal, Input_Segment, Input_Bulk, DropCrit, DropCritNormalStressMultiplier)


    ReturnStiffnessMatrixShear=StiffnessMatrixShearOriginal
    ReturnStiffnessMatrixNormal=StiffnessMatrixNormalOriginal
    ReturnInput_Segment=Input_Segment
    FaultCount1=size(Input_Segment,1)
    LoadingFaultCount=sum(Input_Segment[:,16].>0)
    StableFaultcount=sum(Input_Segment[:,9] - Input_Segment[:,10]  .>0)
    #StableFaultcount=0;#sum(Input_Segment[:,9] - Input_Segment[:,10]  .>0)
    StableFaultcount=sum(Input_Segment[:,9] - Input_Segment[:,10]  .>0)
    #println(StableFaultcount)
    # println(StableFaultcount)
    FaultLengthStrike=Input_Bulk[:,4]
    FaultSegmentToBulk=round.(Int,Input_Segment[:,17])
    DropijPair=zeros(Int,2,1)
    DropFault=0
    DropCount=0
    for i=1:FaultCount1-LoadingFaultCount - StableFaultcount
          for j=1:FaultCount1-LoadingFaultCount - StableFaultcount
            #if StiffnessMatrixShearOriginal[i,j]*StiffnessMatrixShearOriginal[j,i]/
            #        StiffnessMatrixShearOriginal[i,i]/StiffnessMatrixShearOriginal[j,j] > DropCrit
            #if StiffnessMatrixShearOriginal[i,j]/ StiffnessMatrixShearOriginal[i,i] > DropCrit


            InterAction=(StiffnessMatrixShearOriginal[i,j] - DropCritNormalStressMultiplier * StiffnessMatrixNormalOriginal[i,j] ) / StiffnessMatrixShearOriginal[i,i] *
                (StiffnessMatrixShearOriginal[j,i] - DropCritNormalStressMultiplier * StiffnessMatrixNormalOriginal[j,i] ) / StiffnessMatrixShearOriginal[j,j]

            #InterAction=zeros(4)
            #InterAction[1]=(StiffnessMatrixShearOriginal[i,j]+StiffnessMatrixShearOriginal[j,i] ) / StiffnessMatrixShearOriginal[i,i]
            #InterAction[2]=(StiffnessMatrixShearOriginal[j,i]+StiffnessMatrixShearOriginal[i,j] ) / StiffnessMatrixShearOriginal[j,j]

            #InterAction[3]=(DropCritNormalStressMultiplier * StiffnessMatrixNormalOriginal[i,j] + DropCritNormalStressMultiplier * StiffnessMatrixNormalOriginal[j,i]) / StiffnessMatrixShearOriginal[i,i]
            #InterAction[4]=(DropCritNormalStressMultiplier * StiffnessMatrixNormalOriginal[j,i] + DropCritNormalStressMultiplier * StiffnessMatrixNormalOriginal[i,j] ) / StiffnessMatrixShearOriginal[j,j]


            if InterAction > DropCrit
                
                println(i," ", j, " ", InterAction)
                DropCount=DropCount+1
                IandJ=[i,j]
                DropijPair=[DropijPair;IandJ]
                
                if FaultLengthStrike[FaultSegmentToBulk[i]]<FaultLengthStrike[FaultSegmentToBulk[j]]
                    DropFault=[DropFault;i]
                else
                    DropFault=[DropFault;j]
                end            
            end
            
            #if StiffnessMatrixNormalOriginal[i,j]*StiffnessMatrixNormalOriginal[j,i]/
            #        StiffnessMatrixShearOriginal[i,i]/StiffnessMatrixShearOriginal[j,j] > DropCritNormal
            #if StiffnessMatrixNormalOriginal[i,j]/ StiffnessMatrixShearOriginal[i,i]> DropCritNormal
            #    DropCount=DropCount+1
            #    if FaultLengthStrike[FaultSegmentToBulk[i]]<FaultLengthStrike[FaultSegmentToBulk[j]]
            #        DropFault=[DropFault;i]
            #    else
            #        DropFault=[DropFault;j]
            #    end            
            #end
    
    
        end
    end

    if DropCount>0
    DropijPair=DropijPair[2:end,:]
    DropFault=DropFault[2:end,:]
    DropFault=unique(DropFault)
    DropFault=sort(DropFault, rev=true)
    end
    
    println("Drop Faults : ", DropFault)
    
    for i in eachindex(DropFault)
        DropIndex=DropFault[i]
    
        ReturnStiffnessMatrixShear=ReturnStiffnessMatrixShear[1:end .!=DropIndex,:]
        ReturnStiffnessMatrixNormal=ReturnStiffnessMatrixNormal[1:end .!=DropIndex,:]
        ReturnInput_Segment=ReturnInput_Segment[1:end .!=DropIndex,:]

    end
    
    for i in eachindex(DropFault)
        DropIndex=DropFault[i]
    
        ReturnStiffnessMatrixShear=ReturnStiffnessMatrixShear[:, 1:end .!=DropIndex]
        ReturnStiffnessMatrixNormal=ReturnStiffnessMatrixNormal[:,1:end .!=DropIndex]
    end
    
    ReturnFaultCount=size(ReturnInput_Segment,1)
    if length(DropFault) > 0 & maximum(DropFault) > 0
        println("Final Drop Count: ", length(DropFault))
    else 
        println("Final Drop Count: 0")
    end
    println("Return Fault Count: ", ReturnFaultCount)

    return ReturnStiffnessMatrixShear, ReturnStiffnessMatrixNormal, ReturnInput_Segment
end




function SaveResults(StiffnessMatrixShearOriginal, StiffnessMatrixNormalOriginal, ReducedInput_Segment, NormalStiffnessZero, 
    OutputFileName, ShearModulus, PoissonRatio, RockDensity, Switch_StrikeSlip_or_ReverseNormal, MinimumNS);



    FaultCenter=ReducedInput_Segment[:,1:3]
    FaultLengthStrike=ReducedInput_Segment[:,4]
    FaultLengthDip=ReducedInput_Segment[:,5]
    FaultStrikeAngle=ReducedInput_Segment[:,6]
    FaultDipAngle=ReducedInput_Segment[:,7]
    FaultLLRR=ReducedInput_Segment[:,8]
    Fault_a=ReducedInput_Segment[:,9]
    Fault_b=ReducedInput_Segment[:,10]
    Fault_Dc=ReducedInput_Segment[:,11]
    Fault_Theta_i=ReducedInput_Segment[:,12]
    Fault_V_i=ReducedInput_Segment[:,13]
    Fault_Friction_i=ReducedInput_Segment[:,14]
    Fault_NormalStress=ReducedInput_Segment[:,15]
    Fault_V_Const=ReducedInput_Segment[:,16]
    Fault_BulkIndex=ReducedInput_Segment[:,17]
    FaultLengthStrike_Bulk=ReducedInput_Segment[:,18]
    FaultLengthDip_Bulk=ReducedInput_Segment[:,19]

    FaultCount=size(ReducedInput_Segment,1)
    LoadingFaultCount=length(Fault_V_Const[Fault_V_Const.>0])
    FaultMass=RockDensity*(FaultLengthStrike_Bulk+FaultLengthDip_Bulk)/2/(1-PoissonRatio)/pi/pi;
    figure(3)
    clf()
    FaultPlot_3D(ReducedInput_Segment[:,1:3],ReducedInput_Segment[:,4], ReducedInput_Segment[:,5], 
        ReducedInput_Segment[:,6], ReducedInput_Segment[:,7], ReducedInput_Segment[:,8])
    
        xlabel("x")
        ylabel("y")

    figure(3).canvas.draw()

    
    ############################### Save Input File ################################
    ######++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++######

    save(OutputFileName, 
    "StiffnessMatrixShear", StiffnessMatrixShearOriginal, "StiffnessMatrixNormal", StiffnessMatrixNormalOriginal, "FaultCenter", FaultCenter,
    "ShearModulus", ShearModulus, "RockDensity", RockDensity, "PoissonRatio", PoissonRatio,
    "FaultLengthStrike", FaultLengthStrike, "FaultLengthDip", FaultLengthDip, "FaultStrikeAngle", FaultStrikeAngle, 
    "FaultDipAngle", FaultDipAngle, "FaultLLRR", FaultLLRR, "Fault_a", Fault_a, "Fault_b", Fault_b, "Fault_Dc", Fault_Dc, 
    "Fault_Theta_i", Fault_Theta_i, "Fault_V_i", Fault_V_i, "Fault_Friction_i", Fault_Friction_i, "Fault_NormalStress", Fault_NormalStress, 
    "Fault_V_Const", Fault_V_Const, "Fault_BulkIndex", Fault_BulkIndex, "FaultLengthStrike_Bulk", FaultLengthStrike_Bulk, 
    "FaultLengthDip_Bulk", FaultLengthDip_Bulk, "FaultCount", FaultCount, "LoadingFaultCount", LoadingFaultCount, "FaultMass", FaultMass,
    "Switch_StrikeSlip_or_ReverseNormal", Switch_StrikeSlip_or_ReverseNormal, "MinimumNormalStress", MinimumNS,
    "NormalStiffnessZero", NormalStiffnessZero)
    println("Saved File Name: ",OutputFileName)




    # open(SegmentedOutputFileName_List, "w") do io
    #     write(io, "Ctr_X\tCtr_Y\tCtr_Z\tSt_L\tDip_L\tStAng\tDipAng\tLR\ta\tb\tDc\tTheta_i\tV_i\tFric_i\tNormSt\tV_Const\tBulkNo\tBulkSL\tBulkDL\n")
    #     writedlm(io, ReducedInput_Segment)
    # end;

    # println("Saved File Name: ",SegmentedOutputFileName_List)

    ########^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^########
    ################################################################################

end    




function BuildMatrixByPartsStrikeSlip(FaultCount, ElementPartRoughCount, Input_Segment,  ShearModulus, PoissonRatio)

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
            Init_S = (i-1)*PartedElementCount+1
            Fin_S = i*PartedElementCount
            Init_R =  (j-1)*PartedElementCount+1
            Fin_R = j*PartedElementCount
            if i == DivisionCount; Fin_S = FaultCount; end
            if j == DivisionCount; Fin_R = FaultCount; end
    
            StiffnessMatrixShearOriginal[Init_R:Fin_R,Init_S:Fin_S], StiffnessMatrixNormalOriginal[Init_R:Fin_R,Init_S:Fin_S] = 
            StiffnessMatrix_ByParts_Calculation_StrikeSlip(Input_Segment[Init_S:Fin_S,:], Input_Segment[Init_R:Fin_R,:], ShearModulus, PoissonRatio,
                                                CurrentPart, TotalParts)
        end
    end
    
    return StiffnessMatrixShearOriginal, StiffnessMatrixNormalOriginal
end


function BuildMatrixByPartsNormalReverse(FaultCount, ElementPartRoughCount, Input_Segment,  ShearModulus, PoissonRatio)

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
            Init_S = (i-1)*PartedElementCount+1
            Fin_S = i*PartedElementCount
            Init_R =  (j-1)*PartedElementCount+1
            Fin_R = j*PartedElementCount
            if i == DivisionCount; Fin_S = FaultCount; end
            if j == DivisionCount; Fin_R = FaultCount; end
    
            StiffnessMatrixShearOriginal[Init_R:Fin_R,Init_S:Fin_S], StiffnessMatrixNormalOriginal[Init_R:Fin_R,Init_S:Fin_S] = 
            StiffnessMatrix_ByParts_Calculation_NormalReverse(Input_Segment[Init_S:Fin_S,:], Input_Segment[Init_R:Fin_R,:], ShearModulus, PoissonRatio,
                                                CurrentPart, TotalParts)
        end
    end
    
    return StiffnessMatrixShearOriginal, StiffnessMatrixNormalOriginal
end