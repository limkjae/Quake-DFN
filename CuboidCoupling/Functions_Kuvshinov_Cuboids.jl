using CSV
using DataFrames

function Get_Cuboid_Vertices(x, y, z, cube_len_x, cube_len_y, cube_len_z)
    vertices = [
        x - cube_len_x/2   y - cube_len_y/2   z - cube_len_z/2;
        x + cube_len_x/2   y - cube_len_y/2   z - cube_len_z/2;
        x - cube_len_x/2   y + cube_len_y/2   z - cube_len_z/2;
        x + cube_len_x/2   y + cube_len_y/2   z - cube_len_z/2;
        x - cube_len_x/2   y - cube_len_y/2   z + cube_len_z/2;
        x + cube_len_x/2   y - cube_len_y/2   z + cube_len_z/2;
        x - cube_len_x/2   y + cube_len_y/2   z + cube_len_z/2;
        x + cube_len_x/2   y + cube_len_y/2   z + cube_len_z/2
    ]
    return vertices
end


function Load_Reservoir_Cuboids(filename)
    df = CSV.File(filename, header=1, delim=',') |> DataFrame
    MatrixData = Matrix(df)

    Cuboids_Count = size(MatrixData, 1)
    Index = MatrixData[:,1]
    Cuboids_Center = MatrixData[:,2:4]
    Cuboids_Lengths = MatrixData[:,5:7]

    return Cuboids_Count, Cuboids_Center, Cuboids_Lengths
end


function Calculate_Cuboids_Vertices(Cuboids_Count, Cuboids_Center, Cuboids_Lengths)
    Cuboids_Vertices = Vector{Matrix{Float64}}(undef, Cuboids_Count)

    for i = 1:Cuboids_Count
        Cuboids_Vertices[i] = Get_Cuboid_Vertices(
            Cuboids_Center[i,1], Cuboids_Center[i,2], Cuboids_Center[i,3],
            Cuboids_Lengths[i,1], Cuboids_Lengths[i,2], Cuboids_Lengths[i,3]
        )
    end

    return Cuboids_Vertices
end


function Calculate_Correspondense_Matrix(block_center_pos, pressure_pos) 
    """
        Calculate_Correspondense_Matrix(block_center_pos, pressure_pos) 
    
        Return correspondense matrix {C} 
    C[i,j] means the ith cube's pressure change caused by the jth pressure point
    """
    tolerance = 500  #dim change to meters

    correspondance_matrix = zeros( (size(block_center_pos)[1], size(pressure_pos)[1]) )

    for ii in range(1, size(block_center_pos)[1], step = 1)
        #For each reservoir's block, we search the closest point of pressure_pos
        disp_to_each_pressure = sqrt.(  (pressure_pos[:,1] .- block_center_pos[ii,1]).^2 + (pressure_pos[:,2] .- block_center_pos[ii,2]).^2  )
        indx = argmin(disp_to_each_pressure) 
        val  = disp_to_each_pressure[indx] # minimum(disp_to_each_pressure)
        if val < tolerance
            correspondance_matrix[ii, indx] = 1.0
        end
    end

    return correspondance_matrix

end


function fFunc(x,y,z,R)
    return z*atan((x*y)/(z*R)) - x*log(abs(R + y)) - y*log(abs(R + x))
end


function Calculate_PoroDispStress_SingleBlock_FullSpace_GF( ObsPoint, block, block_Cm, ShearModulus  )
    """
    Convention: compression positive + Z downward positive (Left-hand system)
    (same as Flow2Quake and (B.Q. Li et al., 2021) )

    ObsPoint: (1,3), the observation point
    block: (8,3), the vertices of the block
    block_Cm: scalar, the compressibility of the block
    ShearModulus: scalar, the shear modulus of the block
    PossionRatio: scalar, the Possion ratio of the block
    """
    
    pi = 3.1415

    Disp = [0.0, 0.0, 0.0]

    sig_11 = 0.0
    sig_12 = 0.0
    sig_13 = 0.0
    sig_22 = 0.0
    sig_23 = 0.0
    sig_33 = 0.0

    Prefactor_disp = block_Cm/(4*pi)
    Prefactor_sig  = block_Cm*ShearModulus/(2*pi)
    
    for jj in range( start = 1, stop = 8 ) # block:(8,3) Mat
        vertx = block[jj,:] # (1,3)
        xbar  = vertx[1] -  ObsPoint[1]
        ybar  = vertx[2] -  ObsPoint[2]
        ζp = vertx[3] + ObsPoint[3]
        ζm = vertx[3] - ObsPoint[3]
        rp = sqrt(   xbar^2 + ybar^2 + ζp^2    )
        rm = sqrt(   xbar^2 + ybar^2 + ζm^2    )

        Sσ = [-1, 1, 1, -1, 1, -1, -1, 1]

        # -- Displacements --
        # +X  = East, +Y = North, +Z = Down (Left-handed)
        # Compression positive
        Disp[1] +=  Prefactor_disp*Sσ[jj]*fFunc(ybar,ζm,xbar,rm) 
                            
        Disp[2] +=  Prefactor_disp*Sσ[jj]*fFunc(xbar,ζm,ybar,rm)
                            
        Disp[3] +=  - Prefactor_disp*Sσ[jj]*fFunc(xbar,ybar,ζm,rm)
                            
        # -- Stresses --
        # +X  = East, +Y = North, +Z = Down (Left-handed)
        # compression positive
        # σXX
        sig_11 += Prefactor_sig*Sσ[jj]*atan( (xbar*rm)/(ybar*ζm) )
                    
        # σYY
        sig_22 += Prefactor_sig*Sσ[jj]*atan( (ybar*rm)/(xbar*ζm) )

        # σZZ
        sig_33 +=  Prefactor_sig*Sσ[jj]*atan( (ζm*rm)/(xbar*ybar) )

        # σXY
        sig_12 += - Prefactor_sig*Sσ[jj]*log(rm - ζm)
                        
        # σXZ 
        sig_13 += Prefactor_sig*Sσ[jj]*log(rm + ybar)
                        
        # σYZ
        sig_23 += Prefactor_sig*Sσ[jj]*log(rm + xbar)
    end

    
    Sigma = [sig_11, sig_22, sig_33, sig_12, sig_13, sig_23]
    
    return Disp, Sigma 
end





