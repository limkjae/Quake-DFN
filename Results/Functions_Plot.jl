


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

    # FaultPlot_3D_Color_General(FaultCenter[1:FaultCount-LoadingFaultCount,:],FaultLengthStrike[1:FaultCount-LoadingFaultCount], FaultLengthDip[1:FaultCount-LoadingFaultCount],
    # FaultStrikeAngle[1:FaultCount-LoadingFaultCount], FaultDipAngle[1:FaultCount-LoadingFaultCount], FaultLLRR[1:FaultCount-LoadingFaultCount],ReMidPressure[:,1:FaultCount-LoadingFaultCount], 
    # PlotRotation, MinMax_Axis, Transparent)

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
    #MaxDisp=maximum(ReMidDisp)
    #ReMidDisp_Plot=copy(ReMidDisp)
    # for FaultIdx in eachindex(FaultLengthStrike)  
    for FaultIdx = 1: length(FaultLengthStrike)  - LoadingFaultCount
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



