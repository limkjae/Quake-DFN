
function main(StiffnessMatrixShear, StiffnessMatrixNormal, 
    ShearModulus, FaultCount, LoadingFaultCount, Mass,
    a, b, Dc, ThetaI, Vini, FrictionI,
    InitialNormalStress, LoadingRate, 
    TotalStep, RecordStep, SwitchV, TimeStepping, SaveResultFileName,RockDensity,
    FaultCenter,FaultLengthStrike, FaultLengthDip, FaultStrikeAngle, FaultDipAngle, FaultLLRR, SaveStep,
    HMatrixCompress, HMatrix_atol_Shear, HMatrix_atol_Normal, HMatrix_eta, TimeStepOnlyBasedOnUnstablePatch, MinimumNormalStress, Alpha_Evo)
    
    ExternalStressExist=0;

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
    InitialShearStress=InitialNormalStress.*FrictionI;
    
    TotalRecord=Int(round(TotalStep/RecordStep));
    Dt=1;
    ConvergenceCrit=1e-9;


    # define constant variables

    K_Self=zeros(FaultCount)
    Omega=zeros(FaultCount)
    ElasticLoadingShearMatrix=zeros(FaultCount,FaultCount)
    # ElasticLoadingNormal=zeros(FaultCount,FaultCount)

    UnstablePatch=[0;]
    for i=1:FaultCount
        K_Self[i]=abs(StiffnessMatrixShear[i,i]);
        Omega[i]=sqrt(K_Self[i]/Mass[i]);
        for j=1:FaultCount  
            if i==j
                ElasticLoadingShearMatrix[i,j]=0.0;
            else
                ElasticLoadingShearMatrix[i,j]=StiffnessMatrixShear[i,j]/StiffnessMatrixShear[i,i];
                # ElasticLoadingNormal[i,j]=StiffnessMatrixNormal[i,j]/StiffnessMatrixShear[i,i];
            end
        end
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
    if minimum(StiffnessMatrixNormal) - maximum(StiffnessMatrixNormal) == 0
        PlanarFault=1
        println("The fault is planar. Stiffness Matrix Normal is all zero")
    end
    


    Far_Load_Disp_Initial=-(StiffnessMatrixShear\InitialShearStress); # initial load point



    if HMatrixCompress == 1
        Point3D = SVector{3,Float64}
        XCenter = YCenter = [Point3D(FaultCenter[i,:]) for i =1: FaultCount]

        comp = PartialACA(;atol=HMatrix_atol_Shear)
        Xclt = Yclt = ClusterTree(XCenter)
        # adm = StrongAdmissibilityStd()
        adm = StrongAdmissibilityStd(;eta=HMatrix_eta)
        TestX=ones(FaultCount)

        # Build Shear Hmatrix
        Original = ElasticLoadingShearMatrix * TestX
        ElasticLoadingShearMatrix = assemble_hmat(ElasticLoadingShearMatrix,Xclt,Yclt;adm,comp)
        
        println(ElasticLoadingShearMatrix)
        
        Approx = ElasticLoadingShearMatrix * TestX
        Error = abs(maximum((Original-Approx) ./ Original))
        println("Max Shear HMatrix approximation error is ", Error)
        if Error > 1e-8
            println("Warning!!! Your HMatrix_atol may be too large")
            println("Consider smaller HMatrix_atol \n")
        end

        # Build Normal Hmatrix
        if PlanarFault ==0
            XCenter = YCenter = [Point3D(FaultCenter[i,:]) for i =1: FaultCount]
            comp = PartialACA(;atol=HMatrix_atol_Normal)
            Xclt = Yclt = ClusterTree(XCenter)
            adm = StrongAdmissibilityStd()

            Original = StiffnessMatrixNormal * TestX
            StiffnessMatrixNormal = assemble_hmat(StiffnessMatrixNormal,Xclt,Yclt;adm,comp)
            println(StiffnessMatrixNormal)
                
            Approx = StiffnessMatrixNormal * TestX
            Error = abs(maximum((Original-Approx) ./ Original))
            println("Max Normal HMatrix approximation error is ", Error)
                if Error > 1e-6
                    println("Warning!!! Your HMatrix_atol may be too large")
                    println("Consider smaller HMatrix_atol \n")
                end
            else
                println("Planar fault: No H-matrix compression for Normal Stress")
        end

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

    Elastic_Load_Disp_Old=ElasticLoadingShearMatrix*Far_Load_Disp_Initial
    # Elastic_Load_Disp_Old=ElasticLoadingShearMatrix_H*Far_Load_Disp_Initial
    Elastic_Load_Disp=zeros(FaultCount)
    Elastic_Load_Vel=ElasticLoadingShearMatrix*(-VOld)    
    # Elastic_Load_Vel=ElasticLoadingShearMatrix_H*(-VOld)
    Elastic_Load_Vel_Old=copy(Elastic_Load_Vel)
    

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

            while Terminate==0
                Iteration=Iteration+1;
                
                Dt_All=ones(FaultCount,1)*Dt;


                ###################################################
                ##########    External Stress Change    ###########
                if ExternalStressExist ==1
                
                    D_EffStress_Normal, D_EffStress_Shear, D_Pressure = InterpolateFromStressChange(TOld+Dt, FaultCount,
                    ExternalStress_Normal, ExternalStress_Shear,ExternalStress_TimeArray, ExternalStress_Pressure)

                end
                               

                
                ##########    Shear and Normal Stress Change    ###########


                if HMatrixCompress == 1
                    BLAS.set_num_threads(1)
                    mul!(Elastic_Load_Disp, ElasticLoadingShearMatrix, Far_Load_Disp_Initial - DispOld, threads=false)
                    if PlanarFault == 0
                        mul!(EffNormalStressMatrixProduct, StiffnessMatrixNormal, DispOld; threads=false)
                        EffNormalStress_i = EffNormalStressMatrixProduct + InitialNormalStress + D_EffStress_Normal
                    else 
                        EffNormalStress_i = InitialNormalStress + D_EffStress_Normal
                    end

                    BLAS.set_num_threads(Sys.CPU_THREADS)
                else 
                    Elastic_Load_Disp=ElasticLoadingShearMatrix * (Far_Load_Disp_Initial-DispOld) 
                    if PlanarFault == 0
                        EffNormalStressMatrixProduct = StiffnessMatrixNormal * DispOld
                        EffNormalStress_i = EffNormalStressMatrixProduct + InitialNormalStress + D_EffStress_Normal
                    else 
                        EffNormalStress_i = InitialNormalStress + D_EffStress_Normal
                    end

                end
                
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
                    "History_V", History_V, "History_Disp", History_Disp, "History_Pressure", History_Pressure,
                    "History_Time", History_Time, "History_Theta", History_Theta, "History_NormalStress", History_NormalStress,
                    # "History_External_Shear", History_External_Shear, "History_External_Normal", History_External_Normal,
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
            Elastic_Load_Vel_Old=copy(Elastic_Load_Vel)

            
        end
    end

        #println(Return_V)
        #println(Return_Friction)
        #println(Return_Disp)
    

end
