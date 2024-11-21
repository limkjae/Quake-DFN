


function lltokm(latlon1,latlon2)
    radius=6371.0;
    lat1=0*pi/180
    lat2=latlon1*pi/180
    lon1=0*pi/180
    lon2=latlon2*pi/180
    deltaLat=lat2-lat1
    deltaLon=lon2-lon1
    a=sin((deltaLat)/2)^2 + cos(lat1)*cos(lat2) * sin(deltaLon/2)^2
    #c=2*atan2(sqrt(a),sqrt(1-a))
    c=2*atan(sqrt(a)/sqrt(1-a))
    x=deltaLon*cos((lat1+lat2)/2)*radius
    y=deltaLat*radius
    return x, y
end


function FaultPlot_3D_Color_General(FaultCenter,FaultLengthStrike, FaultLengthDip, FaultStrikeAngle,
    FaultDipAngle, FaultLLRR, InputProperty, PlotRotation, MinMax_Axis, ColorMinMax, Transparent, Edge, LoadingFaultCount)

    cm = get_cmap(:jet)    
    buffer=300
    xmin=minimum(FaultCenter[1:end-LoadingFaultCount,1])-buffer
    xmax=maximum(FaultCenter[1:end-LoadingFaultCount,1])+buffer
    ymin=minimum(FaultCenter[1:end-LoadingFaultCount,2])-buffer
    ymax=maximum(FaultCenter[1:end-LoadingFaultCount,2])+buffer
    zmax=0
    zmin=-maximum(FaultCenter[1:end-LoadingFaultCount,3])-buffer

    xCenter=(xmin+xmax)/2
    yCenter=(ymin+ymax)/2
    zCenter=(zmin+zmax)/2
    
    maxlength=maximum([xmax-xmin,ymax-ymin, -zmin+zmax])

    if MinMax_Axis == 0
        xmin=xCenter-maxlength/2
        xmax=xCenter+maxlength/2
        ymin=yCenter-maxlength/2
        ymax=yCenter+maxlength/2
        zmin=zCenter-maxlength/2

    else
        xmin=MinMax_Axis[1,1]
        xmax=MinMax_Axis[1,2]
        ymin=MinMax_Axis[2,1]
        ymax=MinMax_Axis[2,2]
        zmin=MinMax_Axis[3,1]
        zmax=MinMax_Axis[3,2]
    end
    if ColorMinMax==0
    MaxValue=maximum(InputProperty)
    MinValue=minimum(InputProperty)
        if MaxValue == MinValue
            ColorRange = MaxValue*0.1
            MaxValue = MaxValue + ColorRange
            MinValue = MinValue - ColorRange
        end
    else
    MaxValue=ColorMinMax[2]
    MinValue=ColorMinMax[1]
    end

    for FaultIdx = 1: length(FaultLengthStrike)  - LoadingFaultCount

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
        
        fig = figure(1)
        PlotValue=(InputProperty[FaultIdx]-MinValue)/(MaxValue-MinValue)
        art3d = PyObject(PyPlot.art3D)

        verts2 = ([tuple(p1...); tuple(p2...); tuple(p3...); tuple(p4...)],)
        p3c = PyObject(art3d.Poly3DCollection(verts2, linewidths=1))

        ax = subplot(projection="3d")
        pycall(ax.add_collection3d, PyAny, p3c)
        xlim(xmin,xmax )
        ylim(ymin,ymax )
        zlim(zmin,zmax )

        if Transparent == 0
            face_color = [cm(PlotValue)[1], cm(PlotValue)[2],cm(PlotValue)[3],1]#PlotValue]#ReMidVel[PlotStep,FaultIdx]]
            edge_color = [0.5, 0.5, 0.5, 1.0]#ReMidVel[PlotStep,FaultIdx]]
        else 
            face_color = [cm(PlotValue)[1], cm(PlotValue)[2],cm(PlotValue)[3],0.3]#PlotValue]#ReMidVel[PlotStep,FaultIdx]]
            edge_color = [0.5, 0.5, 0.5, 1.0]#ReMidVel[PlotStep,FaultIdx]]

        end
        if Edge == 0
        edge_color = [0 0 0 0]
        end
        
        pycall(p3c.set_facecolor, PyAny, face_color)
        pycall(p3c.set_edgecolor, PyAny, edge_color)
        ax.view_init(PlotRotation[1],PlotRotation[2])

    end
    return MaxValue, MinValue
end




function FaultPlot_3D_Color_General_hsv(FaultCenter,FaultLengthStrike, FaultLengthDip, FaultStrikeAngle,
    FaultDipAngle, FaultLLRR, InputProperty, PlotRotation, MinMax_Axis, ColorMinMax, Transparent, Edge, LoadingFaultCount)

    cm = get_cmap(:hsv)    
    buffer=300
    xmin=minimum(FaultCenter[1:end-LoadingFaultCount,1])-buffer
    xmax=maximum(FaultCenter[1:end-LoadingFaultCount,1])+buffer
    ymin=minimum(FaultCenter[1:end-LoadingFaultCount,2])-buffer
    ymax=maximum(FaultCenter[1:end-LoadingFaultCount,2])+buffer
    zmax=0
    zmin=-maximum(FaultCenter[1:end-LoadingFaultCount,3])-buffer

    xCenter=(xmin+xmax)/2
    yCenter=(ymin+ymax)/2
    zCenter=(zmin+zmax)/2
    
    maxlength=maximum([xmax-xmin,ymax-ymin, -zmin+zmax])

    if MinMax_Axis == 0
        xmin=xCenter-maxlength/2
        xmax=xCenter+maxlength/2
        ymin=yCenter-maxlength/2
        ymax=yCenter+maxlength/2
        zmin=zCenter-maxlength/2

    else
        xmin=MinMax_Axis[1,1]
        xmax=MinMax_Axis[1,2]
        ymin=MinMax_Axis[2,1]
        ymax=MinMax_Axis[2,2]
        zmin=MinMax_Axis[3,1]
        zmax=MinMax_Axis[3,2]
    end
    if ColorMinMax==0
    MaxValue=maximum(InputProperty)
    MinValue=minimum(InputProperty)
        if MaxValue == MinValue
            ColorRange = MaxValue*0.1
            MaxValue = MaxValue + ColorRange
            MinValue = MinValue - ColorRange
        end
    else
    MaxValue=ColorMinMax[2]
    MinValue=ColorMinMax[1]
    end

    for FaultIdx = 1: length(FaultLengthStrike)  - LoadingFaultCount

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
        
        fig = figure(8)
        PlotValue=(InputProperty[FaultIdx]-MinValue)/(MaxValue-MinValue)
        art3d = PyObject(PyPlot.art3D)

        verts2 = ([tuple(p1...); tuple(p2...); tuple(p3...); tuple(p4...)],)
        p3c = PyObject(art3d.Poly3DCollection(verts2, linewidths=1))

        ax = subplot(projection="3d")
        pycall(ax.add_collection3d, PyAny, p3c)
        xlim(xmin,xmax )
        ylim(ymin,ymax )
        zlim(zmin,zmax )

        if Transparent == 0
            face_color = [cm(PlotValue)[1], cm(PlotValue)[2],cm(PlotValue)[3],1]#PlotValue]#ReMidVel[PlotStep,FaultIdx]]
            edge_color = [0.5, 0.5, 0.5, 1.0]#ReMidVel[PlotStep,FaultIdx]]
        else 
            face_color = [cm(PlotValue)[1], cm(PlotValue)[2],cm(PlotValue)[3],0.3]#PlotValue]#ReMidVel[PlotStep,FaultIdx]]
            edge_color = [0.5, 0.5, 0.5, 1.0]#ReMidVel[PlotStep,FaultIdx]]

        end
        if Edge == 0
        edge_color = [0 0 0 0]
        end
        
        pycall(p3c.set_facecolor, PyAny, face_color)
        pycall(p3c.set_edgecolor, PyAny, edge_color)
        ax.view_init(PlotRotation[1],PlotRotation[2])

    end
    return MaxValue, MinValue
end


function FaultPlot_3D_Color_SelectedElements(FaultCenter,FaultLengthStrike, FaultLengthDip, FaultStrikeAngle,
    FaultDipAngle, FaultLLRR, InputProperty, PlotRotation, MinMax_Axis, ColorMinMax, Transparent, SelectedElements)

    # FaultPlot_3D_Color_General(FaultCenter[1:FaultCount-LoadingFaultCount,:],FaultLengthStrike[1:FaultCount-LoadingFaultCount], FaultLengthDip[1:FaultCount-LoadingFaultCount],
    # FaultStrikeAngle[1:FaultCount-LoadingFaultCount], FaultDipAngle[1:FaultCount-LoadingFaultCount], FaultLLRR[1:FaultCount-LoadingFaultCount],ReMidPressure[:,1:FaultCount-LoadingFaultCount], 
    # PlotRotation, MinMax_Axis, Transparent)

    cm = get_cmap(:jet)
    
    buffer=300
    xmin=minimum(FaultCenter[:,1])-buffer
    xmax=maximum(FaultCenter[:,1])+buffer
    ymin=minimum(FaultCenter[:,2])-buffer
    ymax=maximum(FaultCenter[:,2])+buffer
    zmax=0
    zmin=-maximum(FaultCenter[:,3])-buffer

    xCenter=(xmin+xmax)/2
    yCenter=(ymin+ymax)/2
    zCenter=(zmin+zmax)/2
    
    maxlength=maximum([xmax-xmin,ymax-ymin, -zmin+zmax])

    if MinMax_Axis == 0
        xmin=xCenter-maxlength/2
        xmax=xCenter+maxlength/2
        ymin=yCenter-maxlength/2
        ymax=yCenter+maxlength/2
        zmin=zCenter-maxlength/2

    else
        xmin=MinMax_Axis[1,1]
        xmax=MinMax_Axis[1,2]
        ymin=MinMax_Axis[2,1]
        ymax=MinMax_Axis[2,2]
        zmin=MinMax_Axis[3,1]
        zmax=MinMax_Axis[3,2]
    end
    if ColorMinMax==0
    MaxValue=maximum(InputProperty)
    MinValue=minimum(InputProperty)
    else
    MaxValue=ColorMinMax[2]
    MinValue=ColorMinMax[1]
    end
    #MaxDisp=maximum(ReMidDisp)
    #ReMidDisp_Plot=copy(ReMidDisp)
    for FaultIdx in eachindex(FaultLengthStrike)
    # for FaultIdx = 1125
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
        
        fig = figure(1)
        PlotValue=(InputProperty[FaultIdx]-MinValue)/(MaxValue-MinValue)
        
        #if FaultIdx==length(FaultLengthStrike)
        #    ion()
        #end
        #ion() # ioff() for turn off drawnow
        art3d = PyObject(PyPlot.art3D)

        verts2 = ([tuple(p1...); tuple(p2...); tuple(p3...); tuple(p4...)],)
        #p3c = PyObject(art3d.Poly3DCollection(verts2, linewidths=1, alpha=ReMidVel[PlotStep,FaultIdx]/MaxVel))
        #p3c = PyObject(art3d.Poly3DCollection(verts2, linewidths=1, alpha=0.05+ReMidVel[PlotStep,FaultIdx]/MaxVel*0.95)) 
        #p3c = PyObject(art3d.Poly3DCollection(verts2, linewidths=1, alpha=0.05+ReMidVelLog[PlotStep,FaultIdx]/MaxVelLog*0.95))
        p3c = PyObject(art3d.Poly3DCollection(verts2, linewidths=1))
        # ReMidVel[PlotStep,FaultIdx]/MaxVel))

        ax = gca(projection="3d")
        pycall(ax.add_collection3d, PyAny, p3c)
        xlim(xmin,xmax )
        ylim(ymin,ymax )
        zlim(zmin,zmax )

        # if Transparent == 0
        #     face_color = [cm(PlotValue)[1], cm(PlotValue)[2],cm(PlotValue)[3],1]#PlotValue]#ReMidVel[PlotStep,FaultIdx]]
        #     edge_color = [0.5, 0.5, 0.5, 0.1]#ReMidVel[PlotStep,FaultIdx]]
        # else 
        #     face_color = [cm(PlotValue)[1], cm(PlotValue)[2],cm(PlotValue)[3],0.3]#PlotValue]#ReMidVel[PlotStep,FaultIdx]]
        #     edge_color = [0.5, 0.5, 0.5, 0.1]#ReMidVel[PlotStep,FaultIdx]]

        # end
        CheckIfIncluded=0
        for PlotFaultIndex in SelectedElements
            if FaultIdx == PlotFaultIndex
                CheckIfIncluded=CheckIfIncluded+1
            end            
        end        
        if CheckIfIncluded > 0
            face_color = [0 0 0 1]#PlotValue]#ReMidVel[PlotStep,FaultIdx]]
            edge_color = [0 0 0 1]#ReMidVel[PlotStep,FaultIdx]]
        else 
            face_color = [0.1 0.1 0.1 0.1]#PlotValue]#ReMidVel[PlotStep,FaultIdx]]
            edge_color = [0.5 0.5 0.5 0.5]#ReMidVel[PlotStep,FaultIdx]]
        end
        pycall(p3c.set_facecolor, PyAny, face_color)
        pycall(p3c.set_edgecolor, PyAny, edge_color)
        ax.view_init(PlotRotation[1],PlotRotation[2])

    end
    return MaxValue, MinValue


end




function get_event_fragments(FileName, FileNameInput)


    ResultTime=load(FileName,"History_Time")
    ResultDisp=load(FileName,"History_Disp")
    ResultV=load(FileName,"History_V")
    FaultCenter =load(FileNameInput,"FaultCenter")
    FaultLengthStrike =load(FileNameInput, "FaultLengthStrike")
    FaultLengthDip =load(FileNameInput, "FaultLengthDip")
    ShearModulus=load(FileNameInput, "ShearModulus")



    # DispRate = zeros(size(ResultDisp,1),size(ResultDisp,2))
    EventCount_Fragment=1
    EventTime_Fragment=[0;]
    EventMoment_Fragment=[0;]
    EventLocation_Fragment=[0 0 0]

    for FaultIdx in eachindex(ResultDisp[1,:])
         EventOnOrOff=0
         for TimeIdx = 2:length(ResultDisp[:,1])
 
              DispRate= (ResultDisp[TimeIdx,FaultIdx] - ResultDisp[TimeIdx-1,FaultIdx] ) / (ResultTime[TimeIdx] - ResultTime[TimeIdx-1])

              if DispRate > DispRateCrits
                   DispIncrement = ResultDisp[TimeIdx,FaultIdx] - ResultDisp[TimeIdx-1,FaultIdx]
                   MomentIncrement = DispIncrement * FaultLengthStrike[FaultIdx] * FaultLengthDip[FaultIdx]  * ShearModulus     
                   if EventOnOrOff == 0 # If this is a new event
                        EventOnOrOff=1
                        EventCount_Fragment=EventCount_Fragment+1
                        EventTime_Fragment=[EventTime_Fragment; ResultTime[TimeIdx]]
                        EventLocation_Fragment=[EventLocation_Fragment; FaultCenter[FaultIdx,:]']
                        EventMoment_Fragment=[EventMoment_Fragment; MomentIncrement]
                   else  # If this is not a new event
                        EventMoment_Fragment[EventCount_Fragment] = EventMoment_Fragment[EventCount_Fragment] + MomentIncrement 
                   end
              else
                   EventOnOrOff = 0 
              end
         end
    end


    EventCount_Fragment=EventCount_Fragment-1
    EventTime_Fragment=EventTime_Fragment[2:end]
    EventLocation_Fragment=EventLocation_Fragment[2:end,:]
    EventMoment_Fragment=EventMoment_Fragment[2:end]


    ######### Rearrange the fragments by time order ##########
    SortOrder=sortperm(EventTime_Fragment)
    EventTime_Fragment=[EventTime_Fragment[SortOrder];]
    EventLocation_Fragment=EventLocation_Fragment[SortOrder , :]
    EventMoment_Fragment=[EventMoment_Fragment[SortOrder];]


    return  EventTime_Fragment, EventLocation_Fragment, EventMoment_Fragment, EventCount_Fragment
end




function merge_fragments(EventTime_Fragment, EventLocation_Fragment, EventMoment_Fragment, EventCount_Fragment, MergeTimeCriteria, MergeDistanceCriteria)

    SegEventBackIdx=1
    EventCount_Merged=0
    EventTime_Merged=[0;]
    EventLocation_Merged=[0 0 0]

    EventNumberMerged_forFragment=zeros(Int, EventCount_Fragment)
    SumMoment=0.0
    CurrentBulkNumber=0
    for EventFragmentIdx=1:EventCount_Fragment

         Terminate = 0 
         SegEventBackIdx = 1
         while Terminate == 0 
              CheckFragmentIdx = EventFragmentIdx - SegEventBackIdx # check fragment until given time 
              if EventFragmentIdx == 1 # initial setting
                   EventCount_Merged = 1
                   EventTime_Merged = EventTime_Fragment[EventFragmentIdx]
                   EventLocation_Merged = EventLocation_Fragment[EventFragmentIdx, :]'
                   EventNumberMerged_forFragment[EventFragmentIdx] = EventCount_Merged
                   Terminate = 1 
              elseif CheckFragmentIdx < 1
                        EventCount_Merged = EventCount_Merged + 1
                        EventTime_Merged = [EventTime_Merged; EventTime_Fragment[EventFragmentIdx]]
                        EventLocation_Merged = [EventLocation_Merged ; EventLocation_Fragment[EventFragmentIdx, :]']
                        EventNumberMerged_forFragment[EventFragmentIdx] = EventCount_Merged
                        Terminate = 1 
              elseif  EventTime_Fragment[EventFragmentIdx] -  EventTime_Fragment[CheckFragmentIdx] <  MergeTimeCriteria
                   # if within the time range
                   if norm(EventLocation_Fragment[EventFragmentIdx,:] - EventLocation_Fragment[CheckFragmentIdx,:]) < MergeDistanceCriteria
                        # if within the distance range -> This is same event to this event (time and distance satisfied)
                        EventNumberMerged_forFragment[EventFragmentIdx] = EventNumberMerged_forFragment[CheckFragmentIdx] 
                        Terminate = 1 
                   else
                        #this event is within the time range but not in distance range. keep checking
                        SegEventBackIdx = SegEventBackIdx + 1
                   end
              else
                   # Now time range exceeded without finding merged event
                   # This is a new event
                   EventCount_Merged = EventCount_Merged + 1
                   EventTime_Merged = [EventTime_Merged; EventTime_Fragment[EventFragmentIdx]]
                   EventLocation_Merged = [EventLocation_Merged ; EventLocation_Fragment[EventFragmentIdx, :]']
                   EventNumberMerged_forFragment[EventFragmentIdx] = EventCount_Merged
                   Terminate = 1 
              end
         end
    end


    # EventTime_Merged = EventTime_Merged[2:end]
    # EventLocation_Merged = EventLocation_Merged[2:end,:]


    println("new event count = ", EventCount_Merged)

    EventMoment_Bulk = zeros(EventCount_Merged)

    for EventFragmentIdx=1:EventCount_Fragment
         MergedIdxofThisFragment = EventNumberMerged_forFragment[EventFragmentIdx]
         EventMoment_Bulk[MergedIdxofThisFragment] = EventMoment_Bulk[MergedIdxofThisFragment] + EventMoment_Fragment[EventFragmentIdx]
    end
    EventMomentMagnitude = (log10.(EventMoment_Bulk) .- 9.05) ./ 1.5

    return EventTime_Merged, EventMomentMagnitude,  EventLocation_Merged
end






###################################################################
####### Build Stiffness Matrix  Stirke Slip Vector By Part  #######
###################################################################

function SurfaceDeformaion_StrikeSlip(FaultCountSource, FaultCenterSource, FaultLengthStrikeSource,
    FaultLengthDipSource, FaultStrikeAngleSource, FaultDipAngleSource, FaultLLRRSource, 
    FaultCenterReceiver, ShearModulus, PoissonRatio) 

    # FaultCountSource=size(Input_SegmentSource,1)
    # FaultCenterSource=Input_SegmentSource[:,1:3]
    # FaultLengthStrikeSource=Input_SegmentSource[:,4]
    # FaultLengthDipSource=Input_SegmentSource[:,5]
    # FaultStrikeAngleSource=Input_SegmentSource[:,6]
    # FaultDipAngleSource=Input_SegmentSource[:,7]
    # FaultLLRRSource=Input_SegmentSource[:,8]

    FaultCountReceiver=length(FaultCenterReceiver)
    # FaultCenterReceiver=Input_SegmentReceiver[:,1:3]
    FaultStrikeAngleReceiver= zeros(FaultCountReceiver)
    FaultDipAngleReceiver = zeros(FaultCountReceiver)
    FaultLLRRReceiver = ones(FaultCountReceiver)
    # println(FaultCountSource, "  ", FaultCountReceiver)

    DisplacementX = zeros(FaultCountReceiver)
    DisplacementY = zeros(FaultCountReceiver)
    DisplacementZ = zeros(FaultCountReceiver)
    
    # println(CurrentPart,"/",TotalParts)
    
        for SourceIndex=1:FaultCountSource;
            
            println(SourceIndex,"  ",SourceIndex,"/",FaultCountSource)

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
            # println(UX)
            # wait()
            DisplacementX += UX
            DisplacementY += UY
            DisplacementZ += UZ

            # StressXX_SourceFrame=(LameFirstParam*(UXX+UYY+UZZ) + 2*ShearModulus*UXX);
            # StressYY_SourceFrame=(LameFirstParam*(UXX+UYY+UZZ) + 2*ShearModulus*UYY);
            # StressZZ_SourceFrame=(LameFirstParam*(UXX+UYY+UZZ) + 2*ShearModulus*UZZ);
            # StressXY_SourceFrame=(UXY + UYX)*ShearModulus;
            # StressXZ_SourceFrame=(UXZ + UZX)*ShearModulus;
            # StressYZ_SourceFrame=(UYZ + UZY)*ShearModulus;     
            


            # for ReceiverIdx = 1:FaultCountReceiver
            #     Stress_SourceFrame=[StressXX_SourceFrame[ReceiverIdx] StressXY_SourceFrame[ReceiverIdx] StressXZ_SourceFrame[ReceiverIdx]
            #     StressXY_SourceFrame[ReceiverIdx] StressYY_SourceFrame[ReceiverIdx] StressYZ_SourceFrame[ReceiverIdx]
            #     StressXZ_SourceFrame[ReceiverIdx] StressYZ_SourceFrame[ReceiverIdx] StressZZ_SourceFrame[ReceiverIdx]];
        

            #     #######################################################################
            #     ##### Rotate Source Frame Stress to Flat Receiver 

            #     RotationMat_FromReceiver_Strike=
            #     [cosd(-RelativeStrkieAngle[ReceiverIdx]) -sind(-RelativeStrkieAngle[ReceiverIdx])  0
            #     sind(-RelativeStrkieAngle[ReceiverIdx]) cosd(-RelativeStrkieAngle[ReceiverIdx]) 0
            #     0  0  1];

            #     RotationMat_FromReceiver_Dip=
            #     [1 0 0
            #     0 cosd(-ReceiverDipAngle[ReceiverIdx]) -sind(-ReceiverDipAngle[ReceiverIdx])
            #     0 sind(-ReceiverDipAngle[ReceiverIdx]) cosd(-ReceiverDipAngle[ReceiverIdx])]
                
            #     RotationMat_FromReceiver_All = RotationMat_FromReceiver_Dip*RotationMat_FromReceiver_Strike
                                
            #     Stress_Receiver = RotationMat_FromReceiver_All*Stress_SourceFrame*RotationMat_FromReceiver_All'

            #     #######################################################################
            #     ##### Read Normal and Shear Stress Change
            #     StiffnessMatrixNormal[ReceiverIdx,SourceIndex] = - Stress_Receiver[3,3]  # compression is negative
            #     StiffnessMatrixShear[ReceiverIdx,SourceIndex] = - ReceiverLLRR[ReceiverIdx] * Stress_Receiver[1,3]  # right latteral is negative

            #     # println(SourceDipAngle,"  ",Z,"  ",DEPTH, " ", StressZZ_SourceFrame, "  ",Stress_Receiver[3,3])
            # end
        
        end
        print("\033c")                  


        # DisplacementX = UX
        # DisplacementY = UY
        # DisplacementZ = UZ

    return DisplacementX,  DisplacementY, DisplacementZ
end
 
