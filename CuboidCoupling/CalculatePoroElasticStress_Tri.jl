##########################################################################################
#####                                                                               ######
#####   One-way cuboid coupling module with Quake-DFN.                              ######
#####                                                                               ######
#####   Three input file required                                                   ######
#####   1. Input_Cuboids.txt        : Define cuboids' geometry                      ######
#####   2. Input_PorePressure.txt   : Time series of Pressure for each cuboid       ######
#####   3. Input_Temperature.txt    : Time series of Temperature for each cuboid    ######
#####   Run CalculatePoroElasticStress.jl after building input file in quake-DFN    ######
#####                                                                               ######
#####   Written by Qian Shi (qshi2@caltech.edu)                                     ######
#####                                                                               ######
##########################################################################################

using DelimitedFiles
using JLD2
using LinearAlgebra
using Printf
using Statistics

using ProgressBars
using HDF5

# for writing to VTK file for Paraview
using WriteVTK
using StaticArrays

# for loading cuboids and pore pressure files
using CSV
using DataFrames


include("./Functions_Kuvshinov_Cuboids.jl")



function Check_FaultPatch_on_Reservoir(SingleFaultCenter, idx_nearest_cuboid, Cuboids_Center, Cuboids_Length)
   if (abs(SingleFaultCenter[1] - Cuboids_Center[idx_nearest_cuboid,:][1]) <= Cuboids_Length[idx_nearest_cuboid,:][1]/2) &&
      (abs(SingleFaultCenter[2] - Cuboids_Center[idx_nearest_cuboid,:][2]) <= Cuboids_Length[idx_nearest_cuboid,:][2]/2) &&
      (abs(SingleFaultCenter[3] - Cuboids_Center[idx_nearest_cuboid,:][3]) <= Cuboids_Length[idx_nearest_cuboid,:][3]/2) 
        if_faultpatch_on_reservoir = true 
   else
        if_faultpatch_on_reservoir = false
   end

   return if_faultpatch_on_reservoir

end

function Calculate_Nearest_CubeIdx_to_Fault(Cuboids_Center, SingleFaultCenter)
    distances = sum((Cuboids_Center .- SingleFaultCenter') .^ 2, dims=2)
    idx_nearest_cuboid =  argmin(distances)[1]
    return idx_nearest_cuboid
end

function Calculate_Rotation_Matrix_for_Fault( FaultDipAngle, FaultStrikeAngle )
    # Rotate the fault center stress to reference frame to read shear and normal
    RotationMat_FromFault_Strike=
    [cosd(-FaultStrikeAngle) -sind(-FaultStrikeAngle)  0
    sind(-FaultStrikeAngle) cosd(-FaultStrikeAngle) 0
    0  0  1];

    RotationMat_FromFault_Dip=
    [1 0 0
    0 cosd(-FaultDipAngle) -sind(-FaultDipAngle)
    0 sind(-FaultDipAngle) cosd(-FaultDipAngle)]
    
    RotationMat_FromFault = RotationMat_FromFault_Dip*RotationMat_FromFault_Strike;
    
    return RotationMat_FromFault
end

function Rotate_Tensor_to_Stress_on_Fault(RotationMat_FromFault, sigEff_all, FaultRakeAngle)    
    Stress_Fault = RotationMat_FromFault * sigEff_all * RotationMat_FromFault'    
    
    D_Stress_Normal = - Stress_Fault[3,3] # tension is positive (negative to make it compression)
    D_Stress_Shear =  cosd(FaultRakeAngle) * Stress_Fault[1,3] + sind(FaultRakeAngle) * Stress_Fault[2,3]

    return D_Stress_Normal, D_Stress_Shear

end



function main(Input_Fault_File, Input_Cuboids_File, Input_PorePressure_File,Input_Temperature_File, Output_ExternalStress_File)
    # Load Fault parameters
    println("---- Loading Fault and Cuboids  ----")
    RorT = load(Input_Fault_File, "RorT")
    FaultCenter= load(Input_Fault_File, "FaultCenter")
    ShearModulus= load(Input_Fault_File, "ShearModulus")
    PoissonRatio= load(Input_Fault_File, "PoissonRatio")
    FaultStrikeAngle= load(Input_Fault_File, "FaultStrikeAngle")
    FaultDipAngle= load(Input_Fault_File, "FaultDipAngle")
    FaultRakeAngle= load(Input_Fault_File, "FaultRakeAngle")
    FaultCount= load(Input_Fault_File, "FaultCount") 
    Switch_StrikeSlip_or_ReverseNormal = load(Input_Fault_File, "Switch_StrikeSlip_or_ReverseNormal")
    
    
    if RorT =="T"
        FaultCenter[:,3] = -FaultCenter[:,3]
        P1, P2, P3 = load(Input_Fault_File, "P1", "P2", "P3")
   
    end

    println("Fault count is:  ", FaultCount)

    # Elastic Properties
    PwaveModulus = 2 * ShearModulus * (1 - PoissonRatio) / (1 - 2 * PoissonRatio)
    Compressibility = 1/PwaveModulus
    BulkModulus = PwaveModulus - 4/3 * ShearModulus
    ThermalExpansivity = 1e-5 # 1/K
    BoitCoefficient = 1.0 


    # Load Cuboids Positions
    Cuboids_Count, Cuboids_Center, Cuboids_Length = Load_Reservoir_Cuboids(Input_Cuboids_File)
    Cuboids_Vertices = Calculate_Cuboids_Vertices(Cuboids_Count, Cuboids_Center, Cuboids_Length)
    println("Reservoir cuboid count is:  ", Cuboids_Count)

    # Load external time array  
    ExternalStress_TimeArray = readdlm(Input_PorePressure_File, ',' )[1,2:end]
    ExternalStress_TimeArray = Float64.(ExternalStress_TimeArray)
    TimeArrayCount = length(ExternalStress_TimeArray)

    # load external pore pressure change
    PorePressureChange = readdlm(Input_PorePressure_File, ',')[2:end,2:end]'
    PorePressureChange = Float64.(PorePressureChange)

    # load external temperature change
    TemperatureChange = readdlm(Input_Temperature_File, ',')[2:end,2:end]'
    TemperatureChange = Float64.(TemperatureChange)

    # Check the size of PorePressureChange and TemperatureChange should be (TimeArrayCount, FaultCount)
    if size(PorePressureChange, 1) != TimeArrayCount || size(PorePressureChange, 2) != Cuboids_Count
        error("Size of PorePressureChange should be (TimeArrayCount, Cuboids_Count)")
    end

    if size(TemperatureChange, 1) != TimeArrayCount || size(TemperatureChange, 2) != Cuboids_Count
        error("Size of TemperatureChange should be (TimeArrayCount, Cuboids_Count)")
    end


    # External Poro-elastic Stress Change on Fault (initial stresses are assigned in QuickParameterChange.jl)
    ExternalStress_Normal_Poro = zeros(TimeArrayCount,FaultCount)
    ExternalStress_Shear_Poro  = zeros(TimeArrayCount,FaultCount)
    PorePressure_Fault     = zeros(TimeArrayCount,FaultCount)

    PoroStress_11 = zeros(TimeArrayCount,FaultCount)
    PoroStress_22 = zeros(TimeArrayCount,FaultCount)
    PoroStress_33 = zeros(TimeArrayCount,FaultCount)
    PoroStress_12 = zeros(TimeArrayCount,FaultCount)
    PoroStress_13 = zeros(TimeArrayCount,FaultCount)
    PoroStress_23 = zeros(TimeArrayCount,FaultCount)

    PoroDisp_1 = zeros(TimeArrayCount,FaultCount)
    PoroDisp_2 = zeros(TimeArrayCount,FaultCount)
    PoroDisp_3 = zeros(TimeArrayCount,FaultCount)

    PoroDisp_GF_patch_cuboid   = zeros(3, FaultCount, Cuboids_Count) # GF == Green Function
    PoroStress_GF_patch_cuboid = zeros(6, FaultCount, Cuboids_Count)

    println("---- Calculating Green Functions of Poro Elastic Stresses on Fault ----")
    println("---- (", FaultCount, " * ", Cuboids_Count, " pairs) ----")

    for i =  tqdm( 1:FaultCount, unit = "fault patch" )
        for j = 1:Cuboids_Count
            PoroDisp_GF_patch_cuboid[:,i,j], PoroStress_GF_patch_cuboid[:,i,j] = Calculate_PoroDispStress_SingleBlock_FullSpace_GF(FaultCenter[i,:], Cuboids_Vertices[j], Compressibility, ShearModulus  )
        end
    end


    # Preparation work, avoid repeating computation
    println("---- Preparing for Rotatating Stress Tensors ----")
    Idx_Nearest_Cube_All = zeros(Int, FaultCount) # Vector{CartesianIndex}(undef, FaultCount) 
    RotationMat_FromFault_All = zeros(FaultCount, 3, 3) 

    for i = tqdm(1:FaultCount, unit="fault")
        Idx_Nearest_Cube_All[i] = Calculate_Nearest_CubeIdx_to_Fault(Cuboids_Center, FaultCenter[i,:])
        RotationMat_FromFault_All[i, :, :] = Calculate_Rotation_Matrix_for_Fault(FaultDipAngle[i], FaultStrikeAngle[i])
    end



    # triangle Mesh Slip Vector
        UnitVector_Normal = zeros(FaultCount,3)
        UnitVector_StrikeSlip = zeros(FaultCount,3)
        UnitVector_DipSlip = zeros(FaultCount,3)
        UnitVector_Slip = zeros(FaultCount,3)
    for ElemIdx = 1:FaultCount
        P1_i = P1[ElemIdx,:] 
        P2_i = P2[ElemIdx,:] 
        P3_i = P3[ElemIdx,:] 

        UnitVector_Normal[ElemIdx,:]= cross(P2_i-P1_i, P3_i-P1_i) / norm(cross(P2_i-P1_i, P3_i-P1_i))
        UnitVector_StrikeSlip[ElemIdx,:] = cross(UnitVector_Normal[ElemIdx,:], [0, 0, 1]) / 
                                            norm(cross(UnitVector_Normal[ElemIdx,:], [0, 0, 1]) )
        UnitVector_DipSlip[ElemIdx,:] = cross(UnitVector_StrikeSlip[ElemIdx,:], UnitVector_Normal[ElemIdx,:]) /
                                            norm(cross(UnitVector_StrikeSlip[ElemIdx,:], UnitVector_Normal[ElemIdx,:]))
        UnitVector_Slip[ElemIdx,:] = UnitVector_StrikeSlip[ElemIdx,:] * cosd(FaultRakeAngle[ElemIdx]) + 
                                        UnitVector_DipSlip[ElemIdx,:] * sind(FaultRakeAngle[ElemIdx])

    end



    # Rotate Poro Elastic Tensors to Stresses on Fault
    println("---- Rotating Poro Elastic Tensors to Stress on Fault ----")
    for (TimeIdx, Time) in tqdm( enumerate(ExternalStress_TimeArray), unit = " timestep" )
        for i =1:FaultCount
            # Convention: compression positive + Z downward positive (Left-hand system)
            EffectivePressureChange = BoitCoefficient * PorePressureChange[TimeIdx, :] + ThermalExpansivity * BulkModulus * TemperatureChange[TimeIdx, :]
            
            PoroDisp_time_patch   = PoroDisp_GF_patch_cuboid[:,i,:]   * EffectivePressureChange
            PoroStress_time_patch = PoroStress_GF_patch_cuboid[:,i,:] * EffectivePressureChange

            PoroDisp_1[TimeIdx,i]   = PoroDisp_time_patch[1] 
            PoroDisp_2[TimeIdx,i]   = PoroDisp_time_patch[2] 
            PoroDisp_3[TimeIdx,i]   = PoroDisp_time_patch[3] 
            
            PoroStress_11[TimeIdx,i] = PoroStress_time_patch[1] 
            PoroStress_22[TimeIdx,i] = PoroStress_time_patch[2] 
            PoroStress_33[TimeIdx,i] = PoroStress_time_patch[3] 
            PoroStress_12[TimeIdx,i] = PoroStress_time_patch[4] 
            PoroStress_13[TimeIdx,i] = PoroStress_time_patch[5] 
            PoroStress_23[TimeIdx,i] = PoroStress_time_patch[6] 
            
            # Check if fault contact the reservoir 
            if_faultpatch_contact_reservoir  = Check_FaultPatch_on_Reservoir(FaultCenter[i,:], Idx_Nearest_Cube_All[i], Cuboids_Center, Cuboids_Length)
            
            if if_faultpatch_contact_reservoir == true
                PorePressure_Fault[TimeIdx,i] = PorePressureChange[TimeIdx, Idx_Nearest_Cube_All[i]]
            else
                PorePressure_Fault[TimeIdx,i] = 0.0
            end

            # Calculate Effective Stress (Change to Tensional stress positive)
            SigEff_all = zeros(3,3)
            SigEff_all = 
            [-PoroStress_11[TimeIdx,i]     -PoroStress_12[TimeIdx,i]   -PoroStress_13[TimeIdx,i]
            -PoroStress_12[TimeIdx,i]     -PoroStress_22[TimeIdx,i]   -PoroStress_23[TimeIdx,i]
            -PoroStress_13[TimeIdx,i]     -PoroStress_23[TimeIdx,i]   -PoroStress_33[TimeIdx,i] ] +
            [PorePressure_Fault[TimeIdx,i]    0.0            0.0
            0.0            PorePressure_Fault[TimeIdx,i]     0.0
            0.0            0.0            PorePressure_Fault[TimeIdx,i]] 


            if RorT == "R"
                ExternalStress_Normal_Poro[TimeIdx,i], ExternalStress_Shear_Poro[TimeIdx,i] = 
                Rotate_Tensor_to_Stress_on_Fault(RotationMat_FromFault_All[i,:,:], SigEff_all, FaultRakeAngle[i])
            else
                TractionVector = SigEff_all * UnitVector_Normal[i,:]
                ExternalStress_Normal_Poro[TimeIdx,i] = -TractionVector' * UnitVector_Normal[i,:]
                ExternalStress_Shear_Poro[TimeIdx,i]  = TractionVector' * UnitVector_Slip[i,:]
                
                
            end
        end
    end

    # Write to JLD2 File for QuakeDFN solver
    if Switch_StrikeSlip_or_ReverseNormal == 1 
        println("---- Stress calculated in strike-slip orientation ----")
    elseif Switch_StrikeSlip_or_ReverseNormal == 2 
        println("---- Stress calculated in Reverse-Normal orientation ----")
    end

    save(Output_ExternalStress_File,
    "ExternalStress_TimeArray", ExternalStress_TimeArray, 
    "ExternalStress_Normal", ExternalStress_Normal_Poro,
    "ExternalStress_Shear", ExternalStress_Shear_Poro,
    "Pressure",  PorePressure_Fault,
    "PorePressureChange", PorePressureChange,
    "PoroStress_11", PoroStress_11,
    "PoroStress_22", PoroStress_22,
    "PoroStress_33", PoroStress_33,
    "PoroStress_12", PoroStress_12,
    "PoroStress_13", PoroStress_13,
    "PoroStress_23", PoroStress_23,
    "PoroDisp_1",   PoroDisp_1,
    "PoroDisp_2",   PoroDisp_2,
    "PoroDisp_3",   PoroDisp_3)
    # println(ExternalStress_Normal_Poro)
    println("---- JLD2 File Saved: ", Output_ExternalStress_File, " ----")
end


# Input Files 
Input_Fault_File="Input_Discretized.jld2" 
Input_Cuboids_File = "CuboidCoupling/Input_Cuboids.txt"
Input_PorePressure_File = "CuboidCoupling/Input_PorePressure.txt"
Input_Temperature_File = "CuboidCoupling/Input_Temperature.txt"
Output_ExternalStress_File = "Input_ExternalStressChange.jld2"

main(Input_Fault_File, Input_Cuboids_File, Input_PorePressure_File, Input_Temperature_File, Output_ExternalStress_File)