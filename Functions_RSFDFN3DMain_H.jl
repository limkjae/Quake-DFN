
function main_H(ShearModulus, FaultCount, LoadingFaultCount, Mass, NormalStiffnessZero,
    a, b, Dc, ThetaI, Vini, FrictionI,
    InitialNormalStress, LoadingRate, 
    TotalStep, RecordStep, SwitchV, TimeStepping, SaveResultFileName,RockDensity,
    FaultCenter,FaultLengthStrike, FaultLengthDip, FaultStrikeAngle, FaultDipAngle, FaultLLRR, SaveStep,
    TimeStepOnlyBasedOnUnstablePatch, MinimumNormalStress, Alpha_Evo,
    Ranks_Shear, Ranks_Normal, ElementRange_SR, NormalStiffness_H, ShearStiffness_H, ThreadCount)  
    
    ExternalStressExist=0;

    
    ################ Optimizing and Initializing Hmatrix ################
    if ThreadCount == 0
        ThreadCount = Threads.nthreads()
    end
    Epsilon_MaxDiffRatio = 1e-7
    MaxRatioAllowed = 1.5
    MaxIteration = 50

    BlockCount = length(Ranks_Shear)
    println("Shear Stiffness")
    Par_ElementDivision_Shear = ParallelOptimization(ShearStiffness_H, ElementRange_SR, 
                                FaultCount, BlockCount, ThreadCount, MaxRatioAllowed, MaxIteration)
    println("Normal Stiffness")
    Par_ElementDivision_Normal = ParallelOptimization(NormalStiffness_H, ElementRange_SR, 
                                FaultCount, BlockCount, ThreadCount, MaxRatioAllowed, MaxIteration)

    LoadingStiffnessH, K_Self= StiffnessTransitionToLoading(ShearStiffness_H, ElementRange_SR, FaultCount)
    
    println("Finding initial slip distribution")
    InitialShearStress = InitialNormalStress .* FrictionI;
    Far_Load_Disp_Initial = zeros(FaultCount)
    if LoadingFaultCount > 0
        Par_ElementDivision_ShearNoLoading = copy(Par_ElementDivision_Shear)
        Par_ElementDivision_ShearNoLoading[end] = Par_ElementDivision_ShearNoLoading[end]-1
        
        for i=1:ThreadCount -1 
            if Par_ElementDivision_ShearNoLoading[end-i] > Par_ElementDivision_ShearNoLoading[end-i+1] 
                Par_ElementDivision_ShearNoLoading[end-i] = Par_ElementDivision_ShearNoLoading[end-i+1] 
            end
        end
        Far_Load_Disp_Initial[1:end-LoadingFaultCount] = 
                SolveAx_b(LoadingStiffnessH[2:end], K_Self[1:end-LoadingFaultCount], InitialShearStress[1:end-LoadingFaultCount],
                         ElementRange_SR[2:end,:], FaultCount - LoadingFaultCount, Par_ElementDivision_ShearNoLoading, ThreadCount, Epsilon_MaxDiffRatio)
    else
        Far_Load_Disp_Initial = SolveAx_b(LoadingStiffnessH, K_Self, InitialShearStress, ElementRange_SR, FaultCount, 
        Par_ElementDivision_Shear, ThreadCount, Epsilon_MaxDiffRatio)
    end

    # Far_Load_Disp_Initial=-(StiffnessMatrixShear\InitialShearStress); # initial load point
    ##############################################################
    
    if  isfile("Input_ExternalStressChange.jld2")
        
        ExternalStress_Normal = load("Input_ExternalStressChange.jld2", "ExternalStress_Normal")
        ExternalStress_Shear = load("Input_ExternalStressChange.jld2", "ExternalStress_Shear")
        ExternalStress_Pressure = load("Input_ExternalStressChange.jld2", "Pressure")
        ExternalStress_TimeArray = load("Input_ExternalStressChange.jld2", "ExternalStress_TimeArray")
        # "ExternalStress_Normal", ExternalStress_Normal, "ExternalStress_Shear", ExternalStress_Shear, 
        # "ExternalStress_TimeArray", ExternalStress_TimeArray)
        if size(ExternalStress_Normal)[2] == FaultCount
            println("External Stress Change Will be Applied")    
            ExternalStressExist=1;
        else
            println("External Stress Change Dosen't match the fault Count. Will not be applied")    
        end
            
    else
        println("No External Stress Change Detected")
    end




    V0=1e-9;
    
    TotalRecord=Int(round(TotalStep/RecordStep));
    Dt=1;
    ConvergenceCrit=1e-9;


    # define constant variables

    Omega = sqrt.(K_Self ./ Mass);

    # K_Self=zeros(FaultCount)
    # Omega=zeros(FaultCount)
    # ElasticLoadingShearMatrix=zeros(FaultCount,FaultCount)
    # ElasticLoadingNormal=zeros(FaultCount,FaultCount)

    UnstablePatch=[0;]
    for i=1:FaultCount
        # K_Self[i]=abs(StiffnessMatrixShear[i,i]);
        # Omega[i]=sqrt(K_Self[i]/Mass[i]);
        # for j=1:FaultCount  
        #     if i==j
        #         ElasticLoadingShearMatrix[i,j]=0.0;
        #     else
        #         ElasticLoadingShearMatrix[i,j]=StiffnessMatrixShear[i,j]/StiffnessMatrixShear[i,i];
        #         # ElasticLoadingNormal[i,j]=StiffnessMatrixNormal[i,j]/StiffnessMatrixShear[i,i];
        #     end
        # end
        if a[i] - b[i] < 0
            UnstablePatch = [UnstablePatch;i]
        end
    end

    if length(UnstablePatch) == 1
        println("No Unstable Patch")
        TimeStepOnlyBasedOnUnstablePatch = 0
    else
        UnstablePatch = UnstablePatch[2:end]
    end

    PlanarFault=0
    if NormalStiffnessZero == 1
        PlanarFault=1
        println("Stiffness Matrix Normal is all zero")
    end
    


    Friction0=FrictionI-a.*log.(Vini./V0)-b.*log.(ThetaI.*V0./Dc); # Initial friction
    Friction=copy(FrictionI);


    # set zeros
    TOld=0.0;
    Step=0;
    T=0.0;
    TOld=copy(T);

    FLAG_GoodToGo=zeros(FaultCount);

    Accel=zeros(FaultCount);
    DtOld=copy(Dt)
    Disp=zeros(FaultCount);
    DispOld=zeros(FaultCount);
    #DispDelta=copy(Disp)

    Instability=zeros(FaultCount,1);
    InstabilityThistime=zeros(FaultCount,1)

    AccelOld=copy(Accel);

    EffNormalStress=zeros(FaultCount,1)
    EffNormalStress_Old=copy(InitialNormalStress);
    EffNormalStress_i=copy(EffNormalStress)
    EffNormalStressMatrixProduct=zeros(FaultCount)
    D_EffStress_Normal = zeros(FaultCount,1)
    D_EffStress_Shear = zeros(FaultCount,1)
    D_Pressure = zeros(FaultCount,1)
    Dt_All = ones(FaultCount)*Dt

    SolverSwitch=zeros(FaultCount)

    #Total_Loading_Disp_Old=zeros(FaultCount,1)
    #Total_Loading_Disp=zeros(FaultCount,1)
    VOld=copy(Vini);
    V=copy(VOld);
    Theta=copy(ThetaI);
    ThetaOld=copy(Theta);
    FrictionOld=copy(Friction);
    SlowOrFast=0;
    Instability=copy(InstabilityThistime);
    Elastic_Load_Disp_Old = copy(Far_Load_Disp_Initial)
    # Elastic_Load_Disp_Old=ElasticLoadingShearMatrix*Far_Load_Disp_Initial
    # Elastic_Load_Disp_Old=ElasticLoadingShearMatrix_H*Far_Load_Disp_Initial
    Elastic_Load_Disp=zeros(FaultCount)
    # Elastic_Load_Vel=ElasticLoadingShearMatrix*(-VOld)    
    # Elastic_Load_Vel=ElasticLoadingShearMatrix_H*(-VOld)
    # Elastic_Load_Vel_Old=copy(Elastic_Load_Vel)
    

    # Define Histories
    History_Time =zeros(TotalRecord,1);
    History_Disp =zeros(TotalRecord,FaultCount);
    History_V =zeros(TotalRecord,FaultCount);
    History_Dt =zeros(TotalRecord,1);
    History_NormalStress =zeros(TotalRecord,FaultCount);
    History_Theta =zeros(TotalRecord,FaultCount);
    History_Pressure =zeros(TotalRecord,FaultCount);
    History_External_Shear = zeros(TotalRecord,FaultCount);
    History_External_Normal = zeros(TotalRecord,FaultCount);
    
    DtRef=1.0
        @time begin
        for i=1:TotalStep


            FLAG_GoodToGo .= 0
            V=copy(VOld);
            if TimeStepOnlyBasedOnUnstablePatch ==1
                Vmax=maximum(V[UnstablePatch]);
            else
                Vmax=maximum(V[1:length(V)-LoadingFaultCount]);
            end
            

            #Dt = DtAdjust(Dt, maximum(V), SwitchV, MinDt, maximum(Instability))
            DtRef = FunctionDtRef(Vmax, TimeStepping, SlowOrFast)
            if DtRef == TimeStepping[1,4]
                Dt=DtRef;
            elseif Dt<DtRef/1.2
                Dt=Dt*1.2;
            else
                Dt=DtRef;
            end    

            SolverSwitch .= 0
            for FaultIdx=1:FaultCount
                if VOld[FaultIdx] > SwitchV
                    SolverSwitch[FaultIdx] = 1;
                else 
                    SolverSwitch[FaultIdx] = 0;
                end
                
            end


            FaultIdx=0;
            Terminate=0;
            Iteration=0;
            SwitchTime=0;


            NetDisp = Far_Load_Disp_Initial - DispOld
            Elastic_Load_Disp = HmatSolver_Pararllel(NetDisp, LoadingStiffnessH, ElementRange_SR, 
                                 Par_ElementDivision_Shear, ThreadCount, zeros(FaultCount, ThreadCount)) ./ -K_Self

            if PlanarFault == 0
                EffNormalStressMatrixProduct = HmatSolver_Pararllel(DispOld, NormalStiffness_H, ElementRange_SR, 
                                Par_ElementDivision_Normal, ThreadCount, zeros(FaultCount, ThreadCount))  
            end

            while Terminate==0
                Iteration=Iteration+1;
                
                Dt_All=ones(FaultCount,1)*Dt;
                

                ###################################################
                ##########    External Stress Change    ###########
                if ExternalStressExist ==1
                
                    D_EffStress_Normal, D_EffStress_Shear, D_Pressure = InterpolateFromStressChange(TOld+Dt, FaultCount,
                    ExternalStress_Normal, ExternalStress_Shear,ExternalStress_TimeArray, ExternalStress_Pressure)

                end
                                     
                EffNormalStress_i = EffNormalStressMatrixProduct + InitialNormalStress + D_EffStress_Normal
                                
                for iii=1:FaultCount
                    if EffNormalStress_i[iii] < MinimumNormalStress
                        EffNormalStress_i[iii] = MinimumNormalStress;
                    end
                end    
                

                Total_Loading_Disp=Far_Load_Disp_Initial + Elastic_Load_Disp           
                

                ############ Solve It for One Step! ############
                V,Friction,Disp,Theta,EffNormalStress,Dt_All,InstabilityThistime, FLAG_GoodToGo =
                SolveOneTimeStep(ConvergenceCrit,DispOld,FrictionOld,
                ThetaOld,VOld, Friction0,a,b, Dc,V0,K_Self, Dt,Omega,
                Mass,ShearModulus, 
                EffNormalStress_i,Total_Loading_Disp, SolverSwitch,
                V,Friction,Disp,Theta,EffNormalStress,Dt_All,Instability,
                LoadingRate, LoadingFaultCount, FaultCount, RockDensity, D_EffStress_Shear, SwitchV, Alpha_Evo, EffNormalStress_Old)

                ############ End of One Step Solver ############
                # println(findall(x->x==0, FLAG_GoodToGo))
                if minimum(FLAG_GoodToGo) ==0
                # if minimum(Dt_All) < maximum(Dt_All)

                    Dt= minimum(Dt_All);
                    # println(findmin(Dt_All))
                    
                    for FaultIdx=1:FaultCount
                        if InstabilityThistime[FaultIdx]==3
                            #println("Switching to 1 fault index ", FaultIdx, ", Dt ", Dt,", VOld ",VOld[FaultIdx])
                            SolverSwitch[FaultIdx]=1;
                            # Dt=RecommendedTimeStep
                        elseif InstabilityThistime[FaultIdx]==4
                            #println("Switching to 0 fault index ", FaultIdx, ", Dt ", Dt,", VOld ",VOld[FaultIdx])
                            SolverSwitch[FaultIdx]=0;
                            
                            Dt=1e-3
                        end
                    end
                    Terminate=0
                else
                    Terminate=1;
                end

                if Iteration>1e2
                    println("Too many Iterations. Dt is:", Dt)
                    readline()
                end
            end

            
            Dt=minimum(Dt_All);
            
            if rem(i,RecordStep)==0
                Step=Step+1;

                History_Time[Step,1]=TOld;
                History_Dt[Step,1]=Dt;
                History_V[Step,:]=V;
                History_Disp[Step,:]=Disp;
                #History_Dt[Step,1]=Dt;
                #History_Friction[Step,:]=Friction;
                History_Theta[Step,:]=Theta;
                History_NormalStress[Step,:]=EffNormalStress;
                History_Pressure[Step,:] = D_Pressure
                History_External_Shear[Step,:] = D_EffStress_Shear;
                History_External_Normal[Step,:] = D_EffStress_Normal;
                #    History_Pressure(Step,:)=Pressure;
                #History_Accel(Step,:)=Accel;
                #History_Far_Load_Disp(Step,:)=Far_Load_Disp;
                #History_ShearStress(Step,:)=StiffnessMatrixShear*(Far_Load_Disp-Disp);

                @printf("%.5f %.5f   %.3e   %.3e    %i \n", i/TotalStep, T/60/60/24, maximum(V[1:FaultCount-LoadingFaultCount]), Dt ,SlowOrFast )

                if rem(i,SaveStep)==0
                    
                    save(SaveResultFileName, 
                    "History_V", History_V, "History_Disp", History_Disp, "History_Theta", History_Theta,
                    "History_Time", History_Time,  "History_NormalStress", History_NormalStress,
                    # "History_External_Shear", History_External_Shear, "History_External_Normal", History_External_Normal,"History_Pressure", History_Pressure,
                    ) 
                    println("Saved Upto Here")
                    #writedlm("geek.txt", History_V')
                end
                
            end
            
            # if maximum(V[1:end-2])>SwitchV
            if maximum(SolverSwitch[1:FaultCount-LoadingFaultCount]) > 0
                SlowOrFast=1;
            else
                SlowOrFast=0;
            end

            
            T=TOld+Dt;
            Accel=(V.-VOld)./Dt;
            AccelOld=copy(Accel);
            DispOld=copy(Disp);
            ThetaOld=copy(Theta);
            VOld=copy(V);
            FrictionOld=copy(Friction);
            EffNormalStress_Old=copy(EffNormalStress);
            TOld=copy(T);
            Instability=copy(InstabilityThistime);
            DtOld=copy(Dt)
            
            Elastic_Load_Disp_Old=copy(Elastic_Load_Disp)
            # Elastic_Load_Vel_Old=copy(Elastic_Load_Vel)

            
        end
    end

        #println(Return_V)
        #println(Return_Friction)
        #println(Return_Disp)
    

end
