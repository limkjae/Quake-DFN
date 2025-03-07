
function EvolutionLaw(ThetaOld, Dt, V, Dc, Alpha_Evo, EffNormalStress, EffNormalStress_Old, b, Fast1_Slow0, EvolutionDR)

    if EvolutionDR == 1 # Dieterich Law 
        if Fast1_Slow0 == 1 
            Const_T = V/Dc + Alpha_Evo * (EffNormalStress - EffNormalStress_Old) / Dt / b / EffNormalStress
            Theta=(ThetaOld - 1/Const_T)*exp(-V/Dc*Dt) + 1/Const_T; #Dieterich Law Semi analytic with Alpha
        else 
            Theta = (ThetaOld+Dt)/(1 + Dt * V/Dc + Alpha_Evo * Dt * (EffNormalStress - EffNormalStress_Old) / Dt / b / EffNormalStress)
        end
    else
        Theta = ThetaOld-V*ThetaOld/Dc*log(V*ThetaOld/Dc)*Dt; #Ruina Law
    end
    return Theta                
end


# function EvolutionLaw_Slow(ThetaOld, Dt, V, Dc, Alpha_Evo, EffNormalStress, EffNormalStress_Old, b)
#     # Theta = (ThetaOld+Dt)/(1 + Dt * V/Dc + Alpha_Evo * Dt * (EffNormalStress - EffNormalStress_Old) / Dt / b / EffNormalStress)
    
#         Theta = ThetaOld-V*ThetaOld/Dc*log(V*ThetaOld/Dc)*Dt; #Ruina Law
#     return Theta                
# end


function Solver_HighV_D(FaultIdx, ConvergenceCrit,DispOld,FrictionOld,ThetaOld,VOld,
    Friction0,a,b,Dc,V0,K_Self, Dt,Omega,
    Mass,ShearModulus,
    EffNormalStress,Total_Loading_Disp,RockDensity, D_EffStress_Shear, SwitchV, Alpha_Evo, EffNormalStress_Old, EvolutionDR)

FLAG_GoodToGo=1
Iteration=0;
InstabilityThistime=0;
maxDiff=10;
ConvergenceCrit_Iter=copy(ConvergenceCrit);
DV=VOld/100;

V=copy(VOld)/1e10;
Disp=copy(DispOld);
Friction=copy(FrictionOld);
Theta=copy(ThetaOld)

if V<0
    V=1e-100
end
# println("High V Solving Dt is ", Dt)

while maxDiff>ConvergenceCrit_Iter
    Iteration=Iteration+1;



    # 
    VTest=copy(V); # Velocity tested in this NR iteration
    Theta = EvolutionLaw(ThetaOld, Dt, V, Dc, Alpha_Evo, EffNormalStress, EffNormalStress_Old, b, 1, EvolutionDR)

    if Theta < 0        
        Dt=Dt/2; V=copy(VOld); Theta=copy(ThetaOld); Friction=copy(FrictionOld); Disp=copy(DispOld);  maxDiff=0; FLAG_GoodToGo=0; # println("Dt Change 1 ", FaultIdx, " V ", V, " VOld ",VOld, " Dt ", Dt)
        break 
    end


    Friction_RadDamp=Friction0+a*log(V/V0)+b*log(Theta*V0./Dc) + ShearModulus/2/sqrt(ShearModulus/RockDensity)/EffNormalStress*V  - D_EffStress_Shear/EffNormalStress; # Initial friction
    F=Total_Loading_Disp-Friction_RadDamp*EffNormalStress/K_Self;
    Disp=(DispOld-F)*cos(Omega*Dt)+(VOld/Omega)*sin(Omega*Dt)+F; 
    V=(Disp-DispOld)/Dt*2-VOld;
    FOriginal=VTest-V;     


    V=VTest+DV
    VTest2 = copy(V)
    Theta = EvolutionLaw(ThetaOld, Dt, V, Dc, Alpha_Evo, EffNormalStress, EffNormalStress_Old, b, 1, EvolutionDR)

    if Theta < 0        
        Dt=Dt/2; V=copy(VOld); Theta=copy(ThetaOld); Friction=copy(FrictionOld); Disp=copy(DispOld);  maxDiff=0; FLAG_GoodToGo=0; # println("Dt Change 1 ", FaultIdx, " V ", V, " VOld ",VOld, " Dt ", Dt)
        break 
    end
    Friction_RadDamp=Friction0+a*log(V/V0)+b*log(Theta*V0./Dc)+ShearModulus/2/sqrt(ShearModulus/RockDensity)/EffNormalStress*V  - D_EffStress_Shear/EffNormalStress; # Initial friction
    F=Total_Loading_Disp-Friction_RadDamp*EffNormalStress/K_Self;
    Disp=(DispOld-F)*cos(Omega*Dt)+(VOld/Omega)*sin(Omega*Dt)+F; 
    V=(Disp-DispOld)/Dt*2-VOld;

    FUpdated=VTest2-V; # We are testing this Newton Rhapson Function. Lets send this to zero
    
    V=VTest-FOriginal*DV/(FUpdated-FOriginal); # Update velocity

    if V<0
        # println("No convergence at high speed", FaultIdx, " ", V, "  ", Dt)
        if  Dt>1e-5
            # println("Dt Reduced 1_", FaultIdx, " VOld ", VOld, " Dt ", Dt)
            Dt=Dt/2; V=copy(VOld); Theta=copy(ThetaOld); Friction=copy(FrictionOld); Disp=copy(DispOld);
            maxDiff=0;
            FLAG_GoodToGo=0
        else
            Dt=Dt/2;
            InstabilityThistime=4;
            println("Dt Reduced 2_", FaultIdx, " VOld ", VOld, " Dt ", Dt)
            #readline()
            V=copy(VOld);
            Theta=copy(ThetaOld);
            Friction=copy(FrictionOld);
            Disp=copy(DispOld);
            maxDiff=0;
            FLAG_GoodToGo=0
        end

    elseif isnan(V)
        # println("No convergence at high speed", FaultIdx, " ", V, "  ", Dt)
        if  Dt>1e-5
            # println("Dt Reduced 3 ", FaultIdx, " V ", V, " VOld ",VOld, " Dt ", Dt)
            #readline()
            #InstabilityThistime=4;
            Dt=Dt/2;
            V=copy(VOld);
            Theta=copy(ThetaOld);
            Friction=copy(FrictionOld);
            Disp=copy(DispOld);
            maxDiff=0;
            FLAG_GoodToGo=0
        else
            Dt=Dt/2;
            println("Dt Reduced 4 ", FaultIdx, " V ", V, " VOld ",VOld, " Dt ", Dt)
            if Dt<1e-10
            readline()
            end
            InstabilityThistime=4;
            V=copy(VOld);
            Theta=copy(ThetaOld);
            Friction=copy(FrictionOld);
            Disp=copy(DispOld);
            maxDiff=0;
            FLAG_GoodToGo=0

        end

    else 
        Disp=DispOld+(VOld+V)/2*Dt; # Update 
        Theta = EvolutionLaw(ThetaOld, Dt, V, Dc, Alpha_Evo, EffNormalStress, EffNormalStress_Old, b, 1, EvolutionDR)
        # Const_T = V/Dc + Alpha_Evo * (EffNormalStress - EffNormalStress_Old) / Dt / b / EffNormalStress
        # Theta=(ThetaOld - 1/Const_T)*exp(-V/Dc*Dt) + 1/Const_T; #Dieterich Law Semi analytic with Alpha

        if Theta < 0        
            Dt=Dt/2; V=copy(VOld); Theta=copy(ThetaOld); Friction=copy(FrictionOld); Disp=copy(DispOld);  maxDiff=0; FLAG_GoodToGo=0; # println("Dt Change 1 ", FaultIdx, " V ", V, " VOld ",VOld, " Dt ", Dt)
            break 
        end
        Friction=Friction0+a*log(V/V0)+b*log(Theta*V0./Dc); # Initial friction
        maxDiff=abs(VTest./V-1);
        InstabilityThistime=0;
    end

    if rem(Iteration,100)==0;
        if Dt<1e-5
            ConvergenceCrit_Iter=ConvergenceCrit_Iter*10; # Only used when convergence is failed
            InstabilityThistime=2;
            if ConvergenceCrit_Iter > 1e-6
                println("!! Convergence crit increaded at highV to",ConvergenceCrit_Iter, " FaultIdx ", FaultIdx, " V ", V, " Dt ", Dt)
            end
        else
            # println("Dt Recuded (No Convergence) ",ConvergenceCrit_Iter, " ", FaultIdx, " ", V, "  ", Dt)
            Dt=Dt/2; V=copy(VOld); Theta=copy(ThetaOld); Friction=copy(FrictionOld); Disp=copy(DispOld); maxDiff=0; FLAG_GoodToGo=0
        end
    end

    #if ConvergenceCrit_Iter>1e-7
    #end

end

# if isnan(V)
#     println("NaN is here!!!! High V",FaultIdx,"Theta is", Theta)
# end



return V, Friction, Disp, Theta, EffNormalStress, Dt, InstabilityThistime, FLAG_GoodToGo
#Return_Friction=Friction;
#Return_Disp=Disp;
#Return_Theta=Theta;
#Return_EffNormalStress=EffNormalStress;
#Return_Dt=Dt;
#Return_InstabilityThistime= InstabilityThistime;
end


function Solver_LowV_D(FaultIdx, ConvergenceCrit,DispOld,FrictionOld,ThetaOld,VOld,
    Friction0,a,b,Dc,V0,K_Self, Dt,Omega,
    Mass,ShearModulus,
    EffNormalStress,Total_Loading_Disp,RockDensity, D_EffStress_Shear, SwitchV, Alpha_Evo, EffNormalStress_Old, EvolutionDR)


    FLAG_GoodToGo=1
    Iteration=0;
    InstabilityThistime=0;
    maxDiff=10;
    ConvergenceCrit_Iter=copy(ConvergenceCrit);


    V=copy(VOld);
    Disp=copy(DispOld);
    Friction=copy(FrictionOld);
    Theta=copy(ThetaOld);
    while maxDiff>ConvergenceCrit_Iter
        Iteration=Iteration+1;
        VTest=copy(V); # Velocity tested in this NR iteration
        DV=VTest/1e10;

        Accel=(V-VOld)/Dt;
        Theta = EvolutionLaw(ThetaOld, Dt, V, Dc, Alpha_Evo, EffNormalStress, EffNormalStress_Old, b, 0, EvolutionDR)

        if Theta < 0        
            if Dt < SwitchV; InstabilityThistime=3; V=copy(VOld); Theta=copy(ThetaOld); Friction=copy(FrictionOld); Disp=copy(DispOld);  FLAG_GoodToGo=0; # println("Switch 1 ", FaultIdx, " V ", V, " VOld ",VOld, " Dt ", Dt)
            elseif VOld == 0; V=copy(VOld); Theta=copy(ThetaOld); Friction=copy(FrictionOld); Disp=copy(DispOld);  FLAG_GoodToGo=1;  # println("No Switch 1")
            else; Dt=Dt/2; V=copy(VOld); Theta=copy(ThetaOld); Friction=copy(FrictionOld); Disp=copy(DispOld);  maxDiff=0; FLAG_GoodToGo=0; # println("Dt Change 1 ", FaultIdx, " V ", V, " VOld ",VOld, " Dt ", Dt)
            end 
            break 

        elseif isnan(Theta)            
            if Dt < SwitchV; InstabilityThistime=3; V=copy(VOld); Theta=copy(ThetaOld); Friction=copy(FrictionOld); Disp=copy(DispOld);  FLAG_GoodToGo=0; # println("Switch 2 ", FaultIdx, " V ", V, " VOld ",VOld, " Dt ", Dt)
            elseif VOld == 0; V=copy(VOld); Theta=copy(ThetaOld); Friction=copy(FrictionOld); Disp=copy(DispOld);  FLAG_GoodToGo=1;  # println("No Switch 2")
            else; Dt=Dt/2; V=copy(VOld); Theta=copy(ThetaOld); Friction=copy(FrictionOld); Disp=copy(DispOld);  maxDiff=0; FLAG_GoodToGo=0; # println("Dt Change 2 ", FaultIdx, " V ", V, " VOld ",VOld, " Dt ", Dt)
            end 
            break 
        end
        
        Disp=DispOld+(VOld+V)/2*Dt
        Friction= 1/EffNormalStress * (K_Self*(Total_Loading_Disp - Disp)  - ShearModulus/2/sqrt(ShearModulus/RockDensity) * V  + D_EffStress_Shear - Mass*Accel);
        Vel = V0 * exp((Friction-Friction0-b*log(V0*Theta/Dc))/a);
        
        FV = VTest - Vel
        
        # Deviated value for NR
        V=V+DV;
        VTest2 = copy(V)

        Accel=(V-VOld)/Dt;
        Theta = EvolutionLaw(ThetaOld, Dt, V, Dc, Alpha_Evo, EffNormalStress, EffNormalStress_Old, b, 0, EvolutionDR)
        if Theta < 0            
            if Dt < SwitchV; InstabilityThistime=3; V=copy(VOld); Theta=copy(ThetaOld); Friction=copy(FrictionOld); Disp=copy(DispOld);  FLAG_GoodToGo=0; # println("Switch 3 ", FaultIdx, " V ", V, " VOld ",VOld, " Dt ", Dt)
            elseif VOld == 0; V=copy(VOld); Theta=copy(ThetaOld); Friction=copy(FrictionOld); Disp=copy(DispOld);  FLAG_GoodToGo=1;  # println("No Switch 3")
            else; Dt=Dt/2; V=copy(VOld); Theta=copy(ThetaOld); Friction=copy(FrictionOld); Disp=copy(DispOld);  maxDiff=0; FLAG_GoodToGo=0; # println("Dt Change 3 ", FaultIdx, " V ", V, " VOld ",VOld, " Dt ", Dt)
            end 
            break 
        elseif isnan(Theta)     
            if Dt < SwitchV; InstabilityThistime=3; V=copy(VOld); Theta=copy(ThetaOld); Friction=copy(FrictionOld); Disp=copy(DispOld);  FLAG_GoodToGo=0; # println("Switch 4 ", FaultIdx, " V ", V, " VOld ",VOld, " Dt ", Dt)
            elseif VOld == 0; V=copy(VOld); Theta=copy(ThetaOld); Friction=copy(FrictionOld); Disp=copy(DispOld);  FLAG_GoodToGo=1;  # println("No Switch 4")
            else; Dt=Dt/2; V=copy(VOld); Theta=copy(ThetaOld); Friction=copy(FrictionOld); Disp=copy(DispOld);  maxDiff=0; FLAG_GoodToGo=0; # println("Dt Change 4 ", FaultIdx, " V ", V, " VOld ",VOld, " Dt ", Dt)
            end 
            break 
        end
        Disp=DispOld+(VOld+V)/2*Dt
        Friction= 1/EffNormalStress * (K_Self*(Total_Loading_Disp - Disp)  - ShearModulus/2/sqrt(ShearModulus/RockDensity) * V  + D_EffStress_Shear - Mass*Accel);
        Vel = V0 * exp((Friction-Friction0-b*log(V0*Theta/Dc))/a);
        
        FV_DeV = VTest2 - Vel

        V = VTest-FV * DV/(FV_DeV-FV)
        
        if EvolutionDR != 1 && V < 0
            V=exp(log(VTest)-FV*log(VTest2/VTest) /(FV_DeV-FV))
            # println(log(VTest2/VTest))
        end


        Accel=(V-VOld)/Dt;


        Disp=DispOld+(VOld+V)/2*Dt; # Update disp
        Theta = EvolutionLaw(ThetaOld, Dt, V, Dc, Alpha_Evo, EffNormalStress, EffNormalStress_Old, b, 0, EvolutionDR)
        if Theta < 0            
            if Dt < SwitchV; InstabilityThistime=3; V=copy(VOld); Theta=copy(ThetaOld); Friction=copy(FrictionOld); Disp=copy(DispOld);  FLAG_GoodToGo=0; # println("Switch 5 ", FaultIdx, " V ", V, " VOld ",VOld, " Dt ", Dt)
            elseif VOld == 0; V=copy(VOld); Theta=copy(ThetaOld); Friction=copy(FrictionOld); Disp=copy(DispOld);  FLAG_GoodToGo=1;  # println("No Switch 5")
            else; Dt=Dt/2; V=copy(VOld); Theta=copy(ThetaOld); Friction=copy(FrictionOld); Disp=copy(DispOld);  maxDiff=0; FLAG_GoodToGo=0; # println("Dt Change 5 ", FaultIdx, " V ", V, " VOld ",VOld, " Dt ", Dt)
            end 
            break 

        elseif isnan(Theta)   
            if VOld == 0; V=copy(VOld); Theta=copy(ThetaOld); Friction=copy(FrictionOld); Disp=copy(DispOld);  FLAG_GoodToGo=1;  #println("No Switch 6")  
            elseif Dt < 1e-3;               
                if VOld<1e-7; 
                    V=0; Theta=copy(ThetaOld) + Dt; Friction=copy(FrictionOld); Disp=copy(DispOld);  maxDiff=0; FLAG_GoodToGo=1;  
                    println("Fault ", FaultIdx, " calculation failed. Consider reducing time step")
                else InstabilityThistime=3; V=copy(VOld); Theta=copy(ThetaOld); Friction=copy(FrictionOld); Disp=copy(DispOld);  FLAG_GoodToGo=0; 
                    # println("Switch 6 ", FaultIdx, " V ", V, " VOld ",VOld, " Dt ", Dt)
                end
                # InstabilityThistime=3; V=copy(VOld); Theta=copy(ThetaOld); Friction=copy(FrictionOld); Disp=copy(DispOld);  FLAG_GoodToGo=0; # println("Switch 6 ", FaultIdx, " V ", V, " VOld ",VOld, " Dt ", Dt)
            
            else; Dt=Dt/2; V=copy(VOld); Theta=copy(ThetaOld); Friction=copy(FrictionOld); Disp=copy(DispOld);  maxDiff=0; FLAG_GoodToGo=0;  
                # println("Dt Change 6 ", FaultIdx, " V ", V, " VOld ",VOld, " Dt ", Dt)
            end 
            break 
        end
        Friction= 1/EffNormalStress * (K_Self*(Total_Loading_Disp - Disp)  - ShearModulus/2/sqrt(ShearModulus/RockDensity) * V - Mass*Accel);
        # InstabilityThistime=0;
        maxDiff=abs(VTest/V-1);




        if rem(Iteration,200)==0;
            # ConvergenceCrit_Iter=ConvergenceCrit_Iter*10
            if Dt < 1e-4; 
                #println("Switched 7 ", FaultIdx, " ", V, "  ", Dt)    
                InstabilityThistime=3; V=copy(VOld); Theta=copy(ThetaOld); Friction=copy(FrictionOld); Disp=copy(DispOld);  FLAG_GoodToGo=0
                # println("SolverConverted ")
                break
            elseif Dt < 1e1
                ConvergenceCrit_Iter=ConvergenceCrit_Iter*2; # Only used when convergence is failed
                # println("Crit Increased to ", ConvergenceCrit_Iter)
                if ConvergenceCrit_Iter > 1e-6
                    println("!! Convergence crit increaded at low V to ",ConvergenceCrit_Iter, " FaultIdx ", FaultIdx, " V ", V, " Dt ", Dt)
                end
            else
                Dt=Dt/2; V=copy(VOld); Theta=copy(ThetaOld); Friction=copy(FrictionOld); Disp=copy(DispOld);  maxDiff=0; FLAG_GoodToGo=0;
                # println("Dt Change 7 ", FaultIdx, " V ", V, " VOld ",VOld, " Dt ", Dt)
            end

            break
        end
        

    end



return V, Friction, Disp, Theta, EffNormalStress, Dt, InstabilityThistime, FLAG_GoodToGo

end






function SolveOneTimeStep(ConvergenceCrit,DispOld,FrictionOld,
    ThetaOld,VOld, Friction0,a,b, Dc,V0,K_Self, Dt,Omega,
    Mass,ShearModulus, 
    EffNormalStress_i,Total_Loading_Disp, SolverSwitch,
    V,Friction,Disp,Theta,EffNormalStress,Dt_All,InstabilityThistime, 
    LoadingRate, LoadingFaultCount, FaultCount,RockDensity, D_EffStress_Shear, SwitchV, Alpha_Evo, EffNormalStress_Old, EvolutionDR)

    FLAG_GoodToGo=zeros(FaultCount)

    DtOld=copy(Dt)
    for FaultIdx=1:FaultCount-LoadingFaultCount
        if SolverSwitch[FaultIdx] == 1
        # if VOld[FaultIdx] > 1e-3
            
            V[FaultIdx],Friction[FaultIdx],Disp[FaultIdx],Theta[FaultIdx],
            EffNormalStress[FaultIdx],Dt_All[FaultIdx],InstabilityThistime[FaultIdx], FLAG_GoodToGo[FaultIdx] = 
            Solver_HighV_D(FaultIdx, ConvergenceCrit,DispOld[FaultIdx],FrictionOld[FaultIdx],
            ThetaOld[FaultIdx],VOld[FaultIdx], Friction0[FaultIdx],a[FaultIdx],b[FaultIdx],
            Dc[FaultIdx],V0,K_Self[FaultIdx], DtOld, Omega[FaultIdx],
            Mass[FaultIdx],ShearModulus,
            EffNormalStress_i[FaultIdx],Total_Loading_Disp[FaultIdx],RockDensity, D_EffStress_Shear[FaultIdx], SwitchV, Alpha_Evo, EffNormalStress_Old[FaultIdx],EvolutionDR)

        else
            
            V[FaultIdx],Friction[FaultIdx],Disp[FaultIdx],Theta[FaultIdx],
            EffNormalStress[FaultIdx],Dt_All[FaultIdx],InstabilityThistime[FaultIdx], FLAG_GoodToGo[FaultIdx] = 
            Solver_LowV_D(FaultIdx, ConvergenceCrit,DispOld[FaultIdx],FrictionOld[FaultIdx],
            ThetaOld[FaultIdx],VOld[FaultIdx], Friction0[FaultIdx],a[FaultIdx],b[FaultIdx],
            Dc[FaultIdx],V0,K_Self[FaultIdx], DtOld, Omega[FaultIdx],
            Mass[FaultIdx],ShearModulus, 
            EffNormalStress_i[FaultIdx],Total_Loading_Disp[FaultIdx],RockDensity, D_EffStress_Shear[FaultIdx], SwitchV, Alpha_Evo, EffNormalStress_Old[FaultIdx],EvolutionDR)
            
            
        end

    end  
    V[FaultCount-LoadingFaultCount+1:FaultCount].=LoadingRate[FaultCount-LoadingFaultCount+1:FaultCount];
    Friction[FaultCount-LoadingFaultCount+1:FaultCount]=FrictionOld[FaultCount-LoadingFaultCount+1:FaultCount];
    Disp[FaultCount-LoadingFaultCount+1:FaultCount]=DispOld[FaultCount-LoadingFaultCount+1:FaultCount]+V[FaultCount-LoadingFaultCount+1:FaultCount]*Dt;
    Theta[FaultCount-LoadingFaultCount+1:FaultCount]=ThetaOld[FaultCount-LoadingFaultCount+1:FaultCount];
    EffNormalStress[FaultCount-LoadingFaultCount+1:FaultCount].=1e10;
    Dt_All[FaultCount-LoadingFaultCount+1:FaultCount].=Dt;
    InstabilityThistime[FaultCount-LoadingFaultCount+1:FaultCount].=0;
    FLAG_GoodToGo[FaultCount-LoadingFaultCount+1:FaultCount] .= 1
    
return V, Friction, Disp, Theta, EffNormalStress, Dt_All, InstabilityThistime, FLAG_GoodToGo
end




function FunctionDtRef(Vmax, TimeStepping, SlowOrFast)

    if SlowOrFast==1;
        Dt_Lower=TimeStepping[1,3]; # Time step [second] last good working: 2
        Dt_Upper=TimeStepping[1,4]; # Time step [second] last good working: 3
        VTreshLower=TimeStepping[2,3];
        VTreshUpper=TimeStepping[2,4];

        if Vmax<VTreshLower
            DtRef=copy(Dt_Lower);
        elseif Vmax < VTreshUpper
            x=(log10(Vmax)-log10(VTreshLower))/(log10(VTreshUpper)-log10(VTreshLower));
            y=3*x^2-2*x^3;
            LogDt=log10(Dt_Lower)+y*(log10(Dt_Upper)-log10(Dt_Lower));
            DtRef=10^LogDt;        
        else
            DtRef=copy(Dt_Upper);
        end
        
    else    
        Dt_Lower=TimeStepping[1,1]; # Time step [second]
        Dt_Upper=TimeStepping[1,2]; # Time step [second]
        VTreshLower=TimeStepping[2,1];
        VTreshUpper=TimeStepping[2,2];
        
        
        if Vmax<VTreshLower
            DtRef=copy(Dt_Lower);
        elseif  Vmax < VTreshUpper
            x=(log10(Vmax)-log10(VTreshLower))/(log10(VTreshUpper)-log10(VTreshLower));
            y=3*x^2-2*x^3;
            LogDt=log10(Dt_Lower)+y*(log10(Dt_Upper)-log10(Dt_Lower));
            DtRef=10^LogDt;        
        else
            Dt_UpperUpper=TimeStepping[1,3]; # Time step [second]
            VTreshUpperUpper=TimeStepping[2,3];
            x=(log10(Vmax)-log10(VTreshUpper))/(log10(VTreshUpperUpper)-log10(VTreshUpper));
            y=3*x^2-2*x^3;
            LogDt=log10(Dt_Upper)+y*(log10(Dt_UpperUpper)-log10(Dt_Upper));
            DtRef=10^LogDt;        
            # DtRef=copy(Dt_Upper);
        end
    end

    return DtRef
end



function FunctionDTPlot(SwitchV, TimeStepping, RecTimeStep)
    SlowOrFast=0

    testDT=zeros(40)
    Vtest=zeros(40)
    for i=1:40

        V=10^(-i*0.25)
        
        if V>SwitchV
            SlowOrFast=1
        else 
            SlowOrFast=0
        end
        Vtest[i]=V
        testDT[i] = FunctionDtRef(V, TimeStepping, SlowOrFast)

    end
    figure(1)
    plot(log10.(Vtest),log10.(testDT))
    plot([minimum(log10.(Vtest)), maximum(log10.(Vtest))], log10.([RecTimeStep,RecTimeStep]), "r:")
    scatter(log10.(TimeStepping[2,:]),log10.(TimeStepping[1,:]))
    xlabel("log_10(Maximum Velocity)")
    ylabel("log_10(dt)")
end






function InterpolateFromStressChange(Time, FaultCount,
    ExternalStress_Normal, ExternalStress_Shear,ExternalStress_TimeArray, ExternalStress_Pressure)


    TimeArrayIdxCount=length(ExternalStress_TimeArray)

    timeIdx=0;
    for i=1:TimeArrayIdxCount
        if Time >= ExternalStress_TimeArray[i]
            timeIdx=i;
        end
    end
    #println(timeIdx)
    timeArrayIdxLower=timeIdx;
    timeArrayIdxUpper=timeIdx+1;
    if timeIdx==0
        timediff=Time
    else
        timediff=Time-ExternalStress_TimeArray[timeIdx]
    end

    LowerValueEffStress_Normal=zeros(FaultCount)
    LowerValueEffStress_Shear=zeros(FaultCount)
    LowerValueEffStress_Pressure=zeros(FaultCount)

    UpperValueEffStress_Normal=zeros(FaultCount)
    UpperValueEffStress_Shear=zeros(FaultCount)
    UpperValueEffStress_Pressure=zeros(FaultCount)
    D_EffStress_Normal=zeros(FaultCount)
    D_EffStress_Shear=zeros(FaultCount)
    D_EffStress_Pressure=zeros(FaultCount)

    if timeIdx==0
        UpperValueEffStress_Normal = ExternalStress_Normal[timeArrayIdxUpper,:];
        UpperValueEffStress_Shear = ExternalStress_Shear[timeArrayIdxUpper,:];
        UpperValueEffStress_Pressure = ExternalStress_Pressure[timeArrayIdxUpper,:];

        D_EffStress_Normal= timediff/ExternalStress_TimeArray[timeArrayIdxUpper] .*
            (UpperValueEffStress_Normal-LowerValueEffStress_Normal) + LowerValueEffStress_Normal;
        D_EffStress_Shear= timediff/ExternalStress_TimeArray[timeArrayIdxUpper] .*
        (UpperValueEffStress_Shear-LowerValueEffStress_Shear) + LowerValueEffStress_Shear;
        D_EffStress_Pressure= timediff/ExternalStress_TimeArray[timeArrayIdxUpper] .*
        (UpperValueEffStress_Pressure-LowerValueEffStress_Pressure) + LowerValueEffStress_Pressure;

    elseif timeIdx==length(ExternalStress_TimeArray)

        LowerValueEffStress_Normal = ExternalStress_Normal[timeArrayIdxLower,:];
        LowerValueEffStress_Shear = ExternalStress_Shear[timeArrayIdxLower,:];
        LowerValueEffStress_Pressure = ExternalStress_Pressure[timeArrayIdxLower,:];

        D_EffStress_Normal=LowerValueEffStress_Normal
        D_EffStress_Shear=LowerValueEffStress_Shear
        D_EffStress_Pressure=LowerValueEffStress_Pressure

    else

        LowerValueEffStress_Normal = ExternalStress_Normal[timeArrayIdxLower,:];
        LowerValueEffStress_Shear = ExternalStress_Shear[timeArrayIdxLower,:];
        LowerValueEffStress_Pressure = ExternalStress_Pressure[timeArrayIdxLower,:];

        UpperValueEffStress_Normal = ExternalStress_Normal[timeArrayIdxUpper,:];
        UpperValueEffStress_Shear = ExternalStress_Shear[timeArrayIdxUpper,:];
        UpperValueEffStress_Pressure = ExternalStress_Pressure[timeArrayIdxUpper,:];
        
        D_EffStress_Normal= timediff / (ExternalStress_TimeArray[timeArrayIdxUpper]-ExternalStress_TimeArray[timeArrayIdxLower]) .*
            (UpperValueEffStress_Normal-LowerValueEffStress_Normal) + LowerValueEffStress_Normal;
        D_EffStress_Shear= timediff / (ExternalStress_TimeArray[timeArrayIdxUpper]-ExternalStress_TimeArray[timeArrayIdxLower]) .*
        (UpperValueEffStress_Shear-LowerValueEffStress_Shear) + LowerValueEffStress_Shear;
        D_EffStress_Pressure= timediff / (ExternalStress_TimeArray[timeArrayIdxUpper]-ExternalStress_TimeArray[timeArrayIdxLower]) .*
        (UpperValueEffStress_Pressure-LowerValueEffStress_Pressure) + LowerValueEffStress_Pressure;

    end


    return D_EffStress_Normal, D_EffStress_Shear, D_EffStress_Pressure
end





function ReduceTooStrongInteraction_Hmat(StrongInteractionCriteriaMultiple, Admissible,
    FaultCount, ElementRange_SR, ShearStiffness_H, NormalStiffness_H)

    LoadingStiffnessH, K_Self= StiffnessTransitionToLoading(ShearStiffness_H, ElementRange_SR, FaultCount)
    StrongInteractionPair = zeros(Int,1,2)

    BlockCount = length(Admissible)
    for Block = 1:BlockCount
        if Admissible[Block] == 0
            SourceInThisBlock = 0
            for SourceAt = ElementRange_SR[Block,1]:ElementRange_SR[Block,2]
                SourceInThisBlock = SourceInThisBlock+1                
                ReceiverInThisBlock = 0
                for ReceiverAt = ElementRange_SR[Block,3]:ElementRange_SR[Block,4]
                    ReceiverInThisBlock = ReceiverInThisBlock + 1
                    if abs(K_Self[ReceiverAt]) * StrongInteractionCriteriaMultiple < abs( LoadingStiffnessH[Block][ReceiverInThisBlock,SourceInThisBlock]) 
                        RedutionRatio = abs(K_Self[ReceiverAt]) * StrongInteractionCriteriaMultiple / abs(LoadingStiffnessH[Block][ReceiverInThisBlock,SourceInThisBlock])
                        ShearStiffness_H[Block][ReceiverInThisBlock,SourceInThisBlock] = ShearStiffness_H[Block][ReceiverInThisBlock,SourceInThisBlock]  * RedutionRatio
                                    # sign(LoadingStiffnessH[Block][ReceiverInThisBlock,SourceInThisBlock]) * abs(K_Self[ReceiverAt]) * StrongInteractionCriteriaMultiple
                        
                        NormalStiffness_H[Block][ReceiverInThisBlock,SourceInThisBlock] = NormalStiffness_H[Block][ReceiverInThisBlock,SourceInThisBlock] * RedutionRatio
                        StrongInteractionPair = [StrongInteractionPair;  [SourceAt ReceiverAt]]
                    end

                end
            end
        end
    end
    StrongInteractionPair = StrongInteractionPair[2:end,:]

    println(length(StrongInteractionPair[:,1]), "/", FaultCount^2," interaction is smoothened")

    return ShearStiffness_H, NormalStiffness_H
end






function ReduceTooStrongInteraction(StrongInteractionCriteriaMultiple, FaultCount, StiffnessMatrixShear, StiffnessMatrixNormal)


    StrongInteractionPair = zeros(Int,1,2)
    for SourceAt = 1:FaultCount
        for ReceiverAt = 1:FaultCount
            if SourceAt != ReceiverAt
                if abs(StiffnessMatrixShear[ReceiverAt, ReceiverAt]) * StrongInteractionCriteriaMultiple < abs( StiffnessMatrixShear[ReceiverAt, SourceAt]) 
                    RedutionRatio = abs(StiffnessMatrixShear[ReceiverAt, ReceiverAt]) * StrongInteractionCriteriaMultiple /  abs( StiffnessMatrixShear[ReceiverAt, SourceAt]) 
                    StiffnessMatrixShear[ReceiverAt, SourceAt]  = StiffnessMatrixShear[ReceiverAt, SourceAt] * RedutionRatio
                    StiffnessMatrixNormal[ReceiverAt, SourceAt]  = StiffnessMatrixNormal[ReceiverAt, SourceAt] * RedutionRatio
                    StrongInteractionPair = [StrongInteractionPair;  [SourceAt ReceiverAt]]
                end

            end
        end
    end
    StrongInteractionPair = StrongInteractionPair[2:end,:]
    println(length(StrongInteractionPair[:,1]), "/", FaultCount^2," interaction is smoothened")

    return StiffnessMatrixShear, StiffnessMatrixNormal
end