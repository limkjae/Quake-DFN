
using GLMakie
using GeometryBasics
using LinearAlgebra
using JLD2

ResultName="Result_T"
FileName="Results/" * ResultName * ".jld2"
FileNameInput="Results/" * ResultName * "_Input.jld2"

ResultTime, ResultV, ResultDisp, Result_NormalStress, History_Theta, =
load(FileName,"History_Time", "History_V", "History_Disp", "History_NormalStress","History_Theta")
ResultV[ResultV.<=0] .= 1e-100

RorT = load(FileNameInput, "RorT")
FaultCount= load(FileNameInput, "FaultCount")
FaultCenter= load(FileNameInput, "FaultCenter")
FaultLengthStrike= load(FileNameInput, "FaultLengthStrike")
FaultLengthDip= load(FileNameInput, "FaultLengthDip")
FaultStrikeAngle= load(FileNameInput, "FaultStrikeAngle")
FaultDipAngle= load(FileNameInput, "FaultDipAngle")
FaultRakeAngle= load(FileNameInput, "FaultRakeAngle")
LoadingFaultCount= load(FileNameInput, "LoadingFaultCount")
Fault_a= load(FileNameInput, "Fault_a")
Fault_b= load(FileNameInput, "Fault_b")
Fault_Dc= load(FileNameInput, "Fault_Dc")
Fault_Theta_i= load(FileNameInput, "Fault_Theta_i")
Fault_V_i= load(FileNameInput, "Fault_V_i")
Fault_Friction_i= load(FileNameInput, "Fault_Friction_i")
Fault_NormalStress= load(FileNameInput, "Fault_NormalStress")
Fault_BulkIndex= load(FileNameInput, "Fault_BulkIndex")
FaultMass= load(FileNameInput, "FaultMass")

if RorT == "T"
    P1, P2, P3 = load(FileNameInput, "P1", "P2", "P3")
end

# Plot settings
PlotStep = 80
PlotRotation = [30, -100]
Transparent = 0
PolygonAlpha = 1 # 0: fully transparent, 1: fully opaque


@assert 1 <= PlotStep <= length(ResultTime)
NPlot = FaultCount - LoadingFaultCount
@assert NPlot > 0

PlotInput=log10.(ResultV[PlotStep,:]); ColorMinMax=[-12, 0] 
# PlotInput= Result_NormalStress[PlotStep,:] ; ColorMinMax=0
# PlotInput=log10.(ResultPressure[PlotStep,:]); ColorMinMax=[3,6]
# PlotInput= Result_NormalStress[PlotStep,:] -  Fault_NormalStress; ColorMinMax=[-1e6,1e6]
# PlotInput=ResultDisp[PlotStep,:]; ColorMinMax= 0

if ColorMinMax == 0
    ColorMinMax = [minimum(PlotInput), maximum(PlotInput)]
end

# Construct a mesh with three independent vertices per triangle.
vertices = Point3f[]
faces = GeometryBasics.TriangleFace{Int}[]
face_values = Float32[]

if RorT == "T"

    for i in 1:NPlot
        p1 = Point3f(P1[i, 1], P1[i, 2], P1[i, 3])
        p2 = Point3f(P2[i, 1], P2[i, 2], P2[i, 3])
        p3 = Point3f(P3[i, 1], P3[i, 2], P3[i, 3])

        # Skip zero-area or invalid triangles
        if !all(isfinite, (p1, p2, p3))
            continue
        end

        area2 = norm(cross(p2 - p1, p3 - p1))
        if !isfinite(area2) || area2 <= 1e-8
            continue
        end

        n = length(vertices)
        push!(vertices, p1, p2, p3)
        push!(faces, GeometryBasics.TriangleFace{Int}(n+1, n+2, n+3))
        push!(face_values, Float32(PlotInput[i]))
    end
    vertex_values = repeat(face_values, inner = 3)

    @assert !isempty(faces) "No valid triangular fault elements to plot."
else

    for i in 1:NPlot
        RotMatStrike = [
            cosd(FaultStrikeAngle[i]) -sind(FaultStrikeAngle[i]) 0
            sind(FaultStrikeAngle[i])  cosd(FaultStrikeAngle[i]) 0
            0                          0                         1
        ]
        RotMatDip = [
            1 0 0
            0 cosd(FaultDipAngle[i]) -sind(FaultDipAngle[i])
            0 sind(FaultDipAngle[i])  cosd(FaultDipAngle[i])
        ]
        center = [FaultCenter[i, 1], FaultCenter[i, 2], -FaultCenter[i, 3]]
        Ls = FaultLengthStrike[i] / 2
        Ld = FaultLengthDip[i] / 2

        P1 = RotMatStrike * RotMatDip * [ Ls, -Ld, 0.0] + center
        P2 = RotMatStrike * RotMatDip * [-Ls, -Ld, 0.0] + center
        P3 = RotMatStrike * RotMatDip * [-Ls,  Ld, 0.0] + center
        P4 = RotMatStrike * RotMatDip * [ Ls,  Ld, 0.0] + center

        p1 = Point3f(P1...)
        p2 = Point3f(P2...)
        p3 = Point3f(P3...)
        p4 = Point3f(P4...)

        n = length(vertices)

        append!(vertices, [p1, p2, p3, p4])

        # Two triangles per quadrilateral
        push!(faces, GeometryBasics.TriangleFace{Int}(n+1, n+2, n+3))
        push!(faces, GeometryBasics.TriangleFace{Int}(n+1, n+3, n+4))

        # One property value per fault
        push!(face_values, Float32(PlotInput[i]))
    end
        vertex_values = repeat(face_values, inner = 4)
end
# Build the figure
fig = Figure(size = (1100, 800))

ax = Axis3(
    fig[1, 1],
    aspect = :data,
    xlabel = "X",
    ylabel = "Y",
    zlabel = "Z",
    # title = "Log10(Velocity)"
)

# One color value for each vertex.
# Every triangle has three independent vertices, so repeat
# its scalar value three times to give it a uniform color.
mesh!(
    ax,
    vertices,
    faces;
    color = vertex_values,
    colormap = :jet,
    colorrange = (ColorMinMax[1], ColorMinMax[2]),
    alpha = PolygonAlpha,
    transparency = true,
    shading = false
)

Colorbar(
    fig[1, 2],
    limits = (ColorMinMax[1], ColorMinMax[2]),
    colormap = :jet,
    # label = "log10(Velocity)",
    lowclip = :blue,
    highclip = :red
)


ax.azimuth[] = deg2rad(PlotRotation[2])
ax.elevation[] = deg2rad(PlotRotation[1])

display(fig)