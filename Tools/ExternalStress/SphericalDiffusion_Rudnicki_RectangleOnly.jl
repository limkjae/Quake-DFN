
using DelimitedFiles
using JLD2
using LinearAlgebra
using Printf
using SpecialFunctions
using PyPlot
using PyCall
pygui(true)


FlowRate=100 # kg/s
PressureOrigin=[0.0, 0.0,-2000]; # Custom Faults
Permeability = 1e-16;
Viscosity = 0.4e-3;
SkemptonCoeff=0.75;
PoissonRatio_Undrained=0.3;
FluidDensity_Ref = 1e3;

OutputFile="Input_ExternalStressChange.jld2"
LoadingInputFileName="Input_Discretized.jld2" 

############################### Load Input Files ###############################
######++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++######
FaultCenter= load(LoadingInputFileName, "FaultCenter")
ShearModulus= load(LoadingInputFileName, "ShearModulus")
RockDensity= load(LoadingInputFileName, "RockDensity")
PoissonRatio= load(LoadingInputFileName, "PoissonRatio")
FaultLengthStrike= load(LoadingInputFileName, "FaultLengthStrike")
FaultLengthDip= load(LoadingInputFileName, "FaultLengthDip")
FaultStrikeAngle= load(LoadingInputFileName, "FaultStrikeAngle")
FaultDipAngle= load(LoadingInputFileName, "FaultDipAngle")
FaultRakeAngle= load(LoadingInputFileName, "FaultRakeAngle")

Fault_BulkIndex= load(LoadingInputFileName, "Fault_BulkIndex")
FaultLengthStrike_Bulk= load(LoadingInputFileName, "FaultLengthStrike_Bulk")
FaultLengthDip_Bulk= load(LoadingInputFileName, "FaultLengthDip_Bulk")
FaultCount= load(LoadingInputFileName, "FaultCount")
LoadingFaultCount= load(LoadingInputFileName, "LoadingFaultCount")
Switch_StrikeSlip_or_ReverseNormal = load(LoadingInputFileName, "Switch_StrikeSlip_or_ReverseNormal")

########^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^########
################################################################################

TimeCount=150
ExternalStress_TimeArray=zeros(TimeCount)
for i=1:TimeCount
        ExternalStress_TimeArray[i]=0.001*1.2^i
end


function CalculateNSChange(t,FaultCenter, FaultDipAngle, FaultStrikeAngle, FaultRakeAngle, mu, nu, FlowRate, PressureOrigin, 
        Permeability, Viscosity, B, nu_u, rho_0)

    q=FlowRate #Flow rate
    k = Permeability
    eta = Viscosity
    lambd=2*mu*nu/(1-2*nu);
    lambd_u=(2*mu*nu_u)/(1-2*nu_u);
    alpha=3*(nu_u-nu)/(B*(1+nu_u)*(1-2*nu));
    c = k / eta * (lambd_u - lambd) * (lambd + 2*mu) / alpha^2 / (lambd_u + 2* mu) ;
    CenterPoint = [FaultCenter[1],FaultCenter[2],-FaultCenter[3]]
    Distance=CenterPoint-PressureOrigin
    r = norm(CenterPoint-PressureOrigin)
    xi = r/sqrt(c*t);
    g=erf(xi/2)-xi/sqrt(pi)*exp(-xi^2/4);

    sig_11= -q*(lambd_u-lambd)*mu/(4*pi*rho_0*c*r*alpha*(lambd_u+2*mu)) *
            (1*(erfc(xi/2)-2*xi^(-2)*g)+(Distance[1]*Distance[1]/r^2)*(erfc(xi/2)+6*(xi^-2)*g));

    sig_22= -q*(lambd_u-lambd)*mu/(4*pi*rho_0*c*r*alpha*(lambd_u+2*mu)) *
            (1*(erfc(xi/2)-2*xi^(-2)*g)+(Distance[2]*Distance[2]/r^2)*(erfc(xi/2)+6*(xi^-2)*g));

    sig_33= -q*(lambd_u-lambd)*mu/(4*pi*rho_0*c*r*alpha*(lambd_u+2*mu)) *
            (1*(erfc(xi/2)-2*xi^(-2)*g)+(Distance[3]*Distance[3]/r^2)*(erfc(xi/2)+6*(xi^-2)*g));

    sig_12= -q*(lambd_u-lambd)*mu/(4*pi*rho_0*c*r*alpha*(lambd_u+2*mu)) *
            ((Distance[1]*Distance[2]/r^2)*(erfc(xi/2)+6*(xi^-2)*g));

    sig_13= -q*(lambd_u-lambd)*mu/(4*pi*rho_0*c*r*alpha*(lambd_u+2*mu)) *
            ((Distance[1]*Distance[3]/r^2)*(erfc(xi/2)+6*(xi^-2)*g));

    sig_23= -q*(lambd_u-lambd)*mu/(4*pi*rho_0*c*r*alpha*(lambd_u+2*mu)) *
            ((Distance[2]*Distance[3]/r^2)*(erfc(xi/2)+6*(xi^-2)*g));

    Pressure=(q*eta)*erfc(xi/2)/(4*pi*rho_0*r*k);


    # Calculate Effective Stress (Tensional stress is positive)
    sigEff_all=[sig_11 sig_12 sig_13
    sig_12 sig_22 sig_23
    sig_13 sig_23 sig_33] + 
    [Pressure 0 0
    0 Pressure 0
    0 0 Pressure]

    # Rotate the fault center stress to reference frame to read shear and normal
    RotationMat_FromFault_Strike=
    [cosd(-FaultStrikeAngle) -sind(-FaultStrikeAngle)  0
    sind(-FaultStrikeAngle) cosd(-FaultStrikeAngle) 0
    0  0  1];
    RotationMat_FromFault_Dip=
    [1 0 0
    0 cosd(-FaultDipAngle) -sind(-FaultDipAngle)
    0 sind(-FaultDipAngle) cosd(-FaultDipAngle)]
    RotationMat_FromFault_All=RotationMat_FromFault_Dip*RotationMat_FromFault_Strike;
    
    Stress_Fault=RotationMat_FromFault_All * sigEff_all * RotationMat_FromFault_All'    


    D_Stress_Normal = - Stress_Fault[3,3] # tension is positive (negative to make it compression)
    D_Stress_Shear =  cosd(FaultRakeAngle) * Stress_Fault[1,3] + sind(FaultRakeAngle) * Stress_Fault[2,3] 


    return D_Stress_Normal, D_Stress_Shear, Pressure

end



TotalPlotFault=FaultCount
TimeArrayCount=length(ExternalStress_TimeArray)
ExternalStress_Normal=zeros(TimeArrayCount,TotalPlotFault)
ExternalStress_Shear=zeros(TimeArrayCount,TotalPlotFault)
Pressure=zeros(TimeArrayCount,TotalPlotFault)

for TimeIdx in eachindex(ExternalStress_TimeArray)
        Time=ExternalStress_TimeArray[TimeIdx]

        println(Time)
        for i =1:TotalPlotFault
                ExternalStress_Normal[TimeIdx,i], ExternalStress_Shear[TimeIdx,i], Pressure[TimeIdx,i] = 
                        CalculateNSChange(Time,FaultCenter[i,:], FaultDipAngle[i], FaultStrikeAngle[i], FaultRakeAngle[i],  ShearModulus, PoissonRatio, FlowRate, PressureOrigin, 
                        Permeability, Viscosity, SkemptonCoeff, PoissonRatio_Undrained, FluidDensity_Ref)
        end
end

figure(5)
PyPlot.plot(ExternalStress_TimeArray/60/60/24,ExternalStress_Normal)
ylabel("Effective Normal Stress Change (Pa)")
xlabel("Day")
figure(6)
PyPlot.plot(ExternalStress_TimeArray/60/60/24,ExternalStress_Shear)
ylabel("Shear Stress Change (Pa)")
xlabel("Day")


println("File Saved: ", OutputFile)
if Switch_StrikeSlip_or_ReverseNormal == 1 
        println("Stress calculated in strike-slip orientation")
elseif Switch_StrikeSlip_or_ReverseNormal == 2 
        println("Stress calculated in Reverse-Normal orientation")
end

save(OutputFile, 
"ExternalStress_Normal", ExternalStress_Normal, "ExternalStress_Shear", ExternalStress_Shear, "Pressure", Pressure,
"ExternalStress_TimeArray", ExternalStress_TimeArray, "FlowRate", FlowRate, "PressureOrigin", PressureOrigin)
    
