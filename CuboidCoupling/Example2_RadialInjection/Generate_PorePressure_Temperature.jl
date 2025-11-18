using DelimitedFiles
using CSV
using DataFrames
using ProgressBars
using SpecialFunctions: expinti

include("../Functions_Kuvshinov_Cuboids.jl")



function Save_time_data_to_txt(filename::String, time_array::AbstractVector, data::AbstractMatrix)
    open(filename, "w") do io
        write(io, "time(sec), ")
        writedlm(io, time_array', ", ")
        for i in 1:size(data, 2)
            write(io, "cuboid$(i), ")
            writedlm(io, data[:, i]', ", ")
        end
    end
end

function Porepressure_Injection(ExternalStress_TimeArray, Cuboids_Center, InjectionOrigin, MyP_star, MyDiffusivity)
    time_count, blocks_count = size(ExternalStress_TimeArray, 1), size(Cuboids_Center, 1)
    NonUniformPorePressureChange = zeros( time_count, blocks_count )
    Distances_to_injection_origin = zeros( blocks_count )
    
    for idxCube = 1:blocks_count
        Cube_center = Cuboids_Center[idxCube,:]
        Distances_to_injection_origin[idxCube] = norm(Cube_center - InjectionOrigin)
    end

    for (idxTime, Time) in tqdm( enumerate(ExternalStress_TimeArray), unit = " timestep"  )
        for idxCube = 1:blocks_count
            disp_to_injection_origin = Distances_to_injection_origin[idxCube]
            NonUniformPorePressureChange[idxTime, idxCube] = - expinti(  - disp_to_injection_origin^2 / (4*MyDiffusivity*Time)  ) * MyP_star
        end
    end
    return NonUniformPorePressureChange
end 


function main(InputCuboidsFile, OutputPorepressureFile, OutputTemperatureFile)
    # Load Reservoir Centers
    Cuboids_Count, Cuboids_Center, Cuboids_Length = Load_Reservoir_Cuboids(InputCuboidsFile)
    println("---- Cuboids count is:  ", Cuboids_Count, " ----")

    # Time for External Stresses
    ExternalStress_TimeArray = collect(0: 1 : 3) .* (10*365*24*3600) # every 10 years, in seconds
    TimeArrayCount = length(ExternalStress_TimeArray)

    # Pore Pressure Change
    InjectionOrigin = [-1000, 0, 3000]
    MyDiffusivity = 1e-2 # m/second  # 1e-2: 10yr ~ 2000m
    MyP_star = 4e7 # Pa

    PorePressureChange = Porepressure_Injection(ExternalStress_TimeArray, Cuboids_Center, InjectionOrigin, MyP_star, MyDiffusivity)

    # Check the size of PorePressureChange and TemperatureChange should be (TimeArrayCount, FaultCount)
    if size(PorePressureChange, 1) != TimeArrayCount || size(PorePressureChange, 2) != Cuboids_Count
        error("Size of PorePressureChange should be (TimeArrayCount, Cuboids_Count)")
    end

    # Temperature Change
    TemperatureChange = zeros(TimeArrayCount, Cuboids_Count) # No temperature change
    
    # Save Pore Pressure Change to a txt file
    Save_time_data_to_txt(OutputPorepressureFile, ExternalStress_TimeArray, PorePressureChange)
    Save_time_data_to_txt(OutputTemperatureFile, ExternalStress_TimeArray, TemperatureChange)

    println("---- TXT OutputFile Saved: ", OutputPorepressureFile, " ----")
    println("---- TXT OutputFile Saved: ", OutputTemperatureFile, " ----")
end


InputCuboidsFile = "CuboidCoupling/Input_Cuboids.txt"
OutputPorepressureFile = "CuboidCoupling/Input_PorePressure.txt"
OutputTemperatureFile = "CuboidCoupling/Input_Temperature.txt"

main(InputCuboidsFile, OutputPorepressureFile, OutputTemperatureFile)

