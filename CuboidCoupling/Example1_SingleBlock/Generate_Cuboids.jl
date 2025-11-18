# Convention: +X = East, +Y = North, +Z = Down (Left-handed)
using Printf
using CSV
using DataFrames

include("../Functions_Kuvshinov_Cuboids.jl")

function save_cuboids_to_vtk(filename::AbstractString, cubes::Vector{<:AbstractMatrix})
    n_cubes   = length(cubes)
    n_points  = 8 * n_cubes
    n_quads   = 6 * n_cubes
    total_ints = n_quads * (4 + 1)  # (count + 4 vertex indices) per quad

    open(filename, "w") do io
        # Header
        println(io, "# vtk DataFile Version 3.0")
        println(io, "Cubes")
        println(io, "ASCII")
        println(io, "DATASET POLYDATA")

        # Points (flip Z like the PyVista call did with -center[2])
        println(io, "POINTS $n_points float")
        for cube in cubes
            for i in 1:8
                x, y, z = cube[i, 1], cube[i, 2], cube[i, 3]
                @printf(io, "%g %g %g\n", x, y, -z)
            end
        end

        # Faces as quads (6 per cube). Indices are 0-based in VTK legacy format.
        println(io, "POLYGONS $n_quads $total_ints")
        # Vertex ordering compatible with get_cube_vertices_general:
        # 0: (-x,-y,-z), 1: (+x,-y,-z), 2: (-x,+y,-z), 3: (+x,+y,-z),
        # 4: (-x,-y,+z), 5: (+x,-y,+z), 6: (-x,+y,+z), 7: (+x,+y,+z)
        quad_faces = (
            (0, 1, 3, 2),  # bottom (z-)
            (4, 5, 7, 6),  # top (z+)
            (0, 1, 5, 4),  # y-
            (2, 3, 7, 6),  # y+
            (0, 2, 6, 4),  # x-
            (1, 3, 7, 5),  # x+
        )
        for ci in 0:(n_cubes - 1)
            base = ci * 8
            for (i1, i2, i3, i4) in quad_faces
                println(io, "4 $(base+i1) $(base+i2) $(base+i3) $(base+i4)")
            end
        end
    end
end

function save_cuboids_to_txt(filename::AbstractString, cuboids::Vector{<:AbstractMatrix})
    open(filename, "w") do io
        println(io, "Index, CenterX, CenterY, CenterZ, LengthX, LengthY, LengthZ")
    end

    for (idx, cube) in enumerate(cuboids)
        centerx = (cube[1, 1] + cube[end, 1]) / 2
        centery = (cube[1, 2] + cube[end, 2]) / 2
        centerz = (cube[1, 3] + cube[end, 3]) / 2

        x_length = abs(cube[2, 1] - cube[1, 1])
        y_length = abs(cube[3, 2] - cube[1, 2])
        z_length = abs(cube[5, 3] - cube[1, 3])

        open(filename, "a") do io
            @printf(io, "%d, %g, %g, %g, %g, %g, %g\n",
                    idx, centerx, centery, centerz, x_length, y_length, z_length)
        end
    end
end

function main(OutputCubesTXTFilename)
    println("---- Generating Single Cuboid Mesh ----")

    # Reservoir parameters
    Center    = (0.0,0.0, 3000.0)
    LengthX   = 4000.0
    LengthY   = 2000.0
    Thickness = 200.0

    cuboids_all = [
        Get_Cuboid_Vertices(Center[1], Center[2], Center[3], LengthX, LengthY, Thickness)
    ]


    # save_cuboids_to_vtk("CuboidCoupling/Input_Cuboids.vtk", cuboids_all)
    save_cuboids_to_txt(OutputCubesTXTFilename, cuboids_all)

    println("---- TXT file Saved: $(OutputCubesTXTFilename)  ----")
end

OutputCubesTXTFilename = "CuboidCoupling/Input_Cuboids.txt"

main(OutputCubesTXTFilename)
