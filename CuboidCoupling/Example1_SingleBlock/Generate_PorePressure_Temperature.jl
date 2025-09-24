using DelimitedFiles
using CSV
using DataFrames


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




function main(InputCuboidsFile, OutputPorepressureFile, OutputTemperatureFile)
    # Load Reservoir Centers
    Cuboids_Count, Cuboids_Center, Cuboids_Length = Load_Reservoir_Cuboids(InputCuboidsFile)
    println("---- Cuboids count is:  ", Cuboids_Count, " ----")

    # Time for External Stresses
    Time_Count =  2
    ExternalStress_TimeArray = zeros(Time_Count)
    ExternalStress_TimeArray[1]  = 0
    ExternalStress_TimeArray[2] = 10*365*24*3600 # 10 years in seconds

    # Pore Pressure Change
    PorePressureChange = zeros(Time_Count, Cuboids_Count)
    PorePressureChange[1,:] .= 0.0
    PorePressureChange[2,:] .= -30e6 # -30 MPa

    TemperatureChange = zeros(Time_Count, Cuboids_Count)
    TemperatureChange[1,:] .= 0.0
    TemperatureChange[2,:] .= 0.0 # No temperature change

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

