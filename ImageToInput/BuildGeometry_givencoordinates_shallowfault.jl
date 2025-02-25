using DelimitedFiles
using LinearAlgebra


function calculation(i, j, k, L, line)
    Ctr_X = 1000 * (i + k) / 2
    Ctr_Y = 1000 * (j + L) / 2
    if line <= 900
        Ctr_Z = 1000 * 5
        StAng = 180 + atan((L - j) / (k - i)) * (180 / π)
        LR = 1
        Dip_L = 1000 * 10
    else
        Ctr_Z = 1000 * 3.5
        StAng = atan((L - j) / (k - i)) * (180 / π)
        LR = -1
        Dip_L = 1000 * 7
    end
    if StAng > 90
        StAng = StAng - 180
    end
    if StAng < 0
        StAng = StAng + 180
    end
    St_L = 1000 * sqrt((i - k)^2 + (j - L)^2)
    DipAng = 90
    a = 0.002
    b = 0.006
    Dc = 1e-4
    Theta_i = 1e10
    V_i = 1e-15
    Fric_i = 0.6
    Sig0 = 2e6
    SigGrad = 0
    V_Const = 0
    MaxLeng = 1000
    return Ctr_X, Ctr_Y, Ctr_Z, St_L, Dip_L, StAng, DipAng, LR, a, b, Dc, Theta_i, V_i, Fric_i, Sig0, SigGrad, V_Const, MaxLeng
end


data = readdlm("ImageToInput/CollateralPaper_sidefault_segmented_endpoints.txt", ',', Float64)


output = [
    "SwitchSS/RN    ShearMod    PoissonRatio    R_Density   Crit_TooClose   TooCloseNormal_Multiplier   Minimum_NS",
    "1    3.2038e10   0.25    2670.0  1.05    0.6     2.0e6",
    "Ctr_X  Ctr_Y   Ctr_Z   St_L    Dip_L   StAng   DipAng  LR  a   b   Dc  Theta_i     V_i     Fric_i  Sig0    SigGrad     V_Const     MaxLeng"
]


for (index, row) in enumerate(eachrow(data))
    push!(output, join(calculation(row..., index), "\t"))
end


open("Input_BulkFaultGeometry.txt", "w") do f
    for line in output
        println(f, line)
    end
end

println("Done")
