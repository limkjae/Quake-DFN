using DelimitedFiles
using Base
using PyPlot
using PyCall
using JLD2
using LowRankApprox
using Clustering
using LinearAlgebra
@pyimport matplotlib.patches as patches
pygui(true)

############################################################################
############## Hmatrix with Distributed discretization #####################
Hmatrix = true
HowManyDistribution = 5
############################################################################


# TriAlreadyCompiled = @isdefined FunctionRead
# if TriAlreadyCompiled == false
#     include("scripts/Functions_TDstressHS.jl")
#     FunctionRead = 1
# end

include("scripts/Functions_BuildInputFile.jl")
include("scripts/Functions_OKADA3D.jl")
include("scripts/Functions_Plot.jl")
include("scripts/Functions_Hmatrix.jl")




function DiscritizeFullMatrix()

    InputBulkFileName="Input_BulkFaultGeometry.txt"
    OutputFileName="Input_Discretized.jld2"
    HMatrixStructureFile = "Input_HmatrixStructure.jld2"



    Input_Segment, LoadingFaultCount, ShearModulus, PoissonRatio, RockDensity, 
    Switch_StrikeSlip_or_ReverseNormal, DropCrit, DropCritNormalStressMultiplier, MinimumNS, RorT, FaultCount =
        ReadBulkInput(InputBulkFileName)


    ######################## Read Segmented Input File #######################
    FaultCenter, Fault_a, Fault_b, Fault_Dc, Fault_Theta_i, Fault_V_i, Fault_Friction_i, Fault_NormalStress, 
    Fault_V_Const, Fault_BulkIndex, FaultLengthStrike, FaultLengthDip, FaultStrikeAngle, 
    FaultDipAngle, FaultRakeAngle, FaultLengthStrike_Bulk, FaultLengthDip_Bulk, NormalStiffnessZero = 
        ReadSegmentInput(Input_Segment, FaultCount, RorT) 


    ############# Flip if needed and get Unit Vectors (triangle Only) ########
    if RorT == "T"
        P1, P2, P3, UnitVector_Normal, UnitVector_StrikeSlip, UnitVector_DipSlip, UnitVector_Slip = 
            RotVerts_UnitVectors(Input_Segment, FaultCount, FaultRakeAngle) 
    end


    if RorT == "R"
        StiffnessMatrixShear, StiffnessMatrixNormal = 
            BuildMatrixByParts_Even_Rec(FaultCount, Input_Segment,  ShearModulus, PoissonRatio)
    elseif RorT == "T"
        StiffnessMatrixShear, StiffnessMatrixNormal = 
            BuildMatrixByParts_Even_Tri(P1, P2, P3, FaultRakeAngle, FaultCenter, UnitVector_Normal, UnitVector_Slip,
                                FaultCount, ShearModulus, PoissonRatio)                            
    end


    save(OutputFileName, 
    "FaultCenter", FaultCenter,
    "ShearModulus", ShearModulus, "RockDensity", RockDensity, "PoissonRatio", PoissonRatio,
    "FaultLengthStrike", FaultLengthStrike, "FaultLengthDip", FaultLengthDip, "FaultStrikeAngle", FaultStrikeAngle, 
    "FaultDipAngle", FaultDipAngle, "FaultRakeAngle", FaultRakeAngle, "Fault_a", Fault_a, "Fault_b", Fault_b, "Fault_Dc", Fault_Dc, 
    "Fault_Theta_i", Fault_Theta_i, "Fault_V_i", Fault_V_i, "Fault_Friction_i", Fault_Friction_i, "Fault_NormalStress", Fault_NormalStress, 
    "Fault_V_Const", Fault_V_Const, "Fault_BulkIndex", Fault_BulkIndex, "FaultLengthStrike_Bulk", FaultLengthStrike_Bulk, 
    "FaultLengthDip_Bulk", FaultLengthDip_Bulk, "FaultCount", FaultCount, "LoadingFaultCount", LoadingFaultCount, 
    "Switch_StrikeSlip_or_ReverseNormal", Switch_StrikeSlip_or_ReverseNormal, "MinimumNormalStress", MinimumNS,
    "NormalStiffnessZero", NormalStiffnessZero, "RorT", RorT,   "Hmatrix", Hmatrix)
    println("Saved File Name: ",OutputFileName)

    file = jldopen(OutputFileName, "a+")
        write(file, "StiffnessMatrixShear", StiffnessMatrixShear) 
        write(file, "StiffnessMatrixNormal", StiffnessMatrixNormal) 
    close(file)

    ######################## Save Vertices for Triangles ##############################
    if RorT == "T"        
        file = jldopen(OutputFileName, "a+")
        write(file, "P1", P1) 
        write(file, "P2", P2) 
        write(file, "P3", P3) 
        close(file)
    end


end






function DistDisc()


    InputBulkFileName="Input_BulkFaultGeometry.txt"
    HMatrixStructureFile = "Input_HmatrixStructure.jld2"


    Block_Ctr_Diam, Block_Range_Level, Input_Segment, LoadingFaultExist, LoadingFaultCount,
    MinimumElementsToCut, ArrangePoint, Admissible, ElementRange_SR, Switch_StrikeSlip_or_ReverseNormal, 
    ShearModulus, PoissonRatio, RockDensity, DropCrit, DropCritNormalStressMultiplier, MinimumNS, RorT =
        load(HMatrixStructureFile, "Block_Ctr_Diam", "Block_Range_Level", "Input_Segment", "LoadingFaultExist", "LoadingFaultCount",
            "MinimumElementsToCut", "ArrangePoint", "Admissible", "ElementRange_SR", "Switch_StrikeSlip_or_ReverseNormal", 
            "ShearModulus", "PoissonRatio", "RockDensity", "DropCrit", "DropCritNormalStressMultiplier", "MinimumNS", "RorT")


    BlockElementCount = (ElementRange_SR[:,2] - ElementRange_SR[:,1] .+ 1) .* (ElementRange_SR[:,4] - ElementRange_SR[:,3] .+ 1)  

    BlockElementIncrementNorm = cumsum(BlockElementCount) / sum(BlockElementCount)
    DivisionCut = zeros(Int, HowManyDistribution)
    DivisionBlockItoF = zeros(Int, HowManyDistribution,2)
    TotalBlock = length(BlockElementIncrementNorm)


    for BlockIndex in eachindex(BlockElementIncrementNorm)
        for DivisionIndex = 1 : HowManyDistribution - 1
            DivisionCrit =  DivisionIndex / HowManyDistribution 
            if BlockElementIncrementNorm[BlockIndex] < DivisionCrit
                DivisionCut[DivisionIndex] = BlockIndex
            
            end
        end
    end
    DivisionCut[end] = length(BlockElementIncrementNorm)

    for index = 1:HowManyDistribution
        if index == 1
            DivisionBlockItoF[1,:] = [1, DivisionCut[index]] 
        else
            DivisionBlockItoF[index,:] = [DivisionCut[index-1]+1, DivisionCut[index]] 

        end
    end




    for DistributeIndex = 1:HowManyDistribution

    
    if Sys.isapple()
        InputBulkFileNamePart = "$(@__DIR__)/Input_BulkFaultGeometry.txt"
        OutputFileNamePart = "$(@__DIR__)/temp_HMatPart_$(DistributeIndex).jld2"
        HMatrixStructureFilePart = "$(@__DIR__)/Input_HmatrixStructure.jld2"
    else
        InputBulkFileNamePart = "Input_BulkFaultGeometry.txt"
        OutputFileNamePart = "temp_HMatPart_$(DistributeIndex).jld2"
        HMatrixStructureFilePart = "Input_HmatrixStructure.jld2"
    end

    if RorT == "T"
        ExtraUsing = "using TriangularDislocation"
    else
        ExtraUsing = ""
    end


    CodeScript = """


        println("Distributed Discretization part $(DistributeIndex) / $(HowManyDistribution)")
        println("Compiling... This may take a few minutes")
        

        using DelimitedFiles
        using Base
        using PyPlot
        using PyCall
        using JLD2
        using LowRankApprox
        using Clustering
        using LinearAlgebra
        $(ExtraUsing)

        include("scripts/Functions_BuildInputFile.jl")
        include("scripts/Functions_OKADA3D.jl")
        include("scripts/Functions_Plot.jl")
        include("scripts/Functions_Hmatrix.jl")
        
        

        function Discritize()

            BlockI = $(DivisionBlockItoF[DistributeIndex,1])
            BlockF = $(DivisionBlockItoF[DistributeIndex,2])

            InputBulkFileName = "$(InputBulkFileNamePart)"
            OutputFileName =  "$(OutputFileNamePart)"
            HMatrixStructureFile = "$(HMatrixStructureFilePart)"

            ##########################       Hmatrix       #############################
            Tolerance = 1e3 #  Hmatrix compression Tolerance. pascal for 1m slip (More approximaion for higher Tolerance). 

            Block_Ctr_Diam, Block_Range_Level, Input_Segment, LoadingFaultExist, LoadingFaultCount,
            MinimumElementsToCut, ArrangePoint, Admissible, ElementRange_SR, Switch_StrikeSlip_or_ReverseNormal, 
            ShearModulus, PoissonRatio, RockDensity, DropCrit, DropCritNormalStressMultiplier, MinimumNS, RorT =
            load(HMatrixStructureFile, "Block_Ctr_Diam", "Block_Range_Level", "Input_Segment", "LoadingFaultExist", "LoadingFaultCount",
                "MinimumElementsToCut", "ArrangePoint", "Admissible", "ElementRange_SR", "Switch_StrikeSlip_or_ReverseNormal", 
                "ShearModulus", "PoissonRatio", "RockDensity", "DropCrit", "DropCritNormalStressMultiplier", "MinimumNS", "RorT")

            Admissible = Admissible[BlockI:BlockF]
            ElementRange_SR = ElementRange_SR[BlockI:BlockF,:]



            FaultCount = length(Input_Segment[:,1])


            ######################## Read Segmented Input File #######################
            FaultCenter, Fault_a, Fault_b, Fault_Dc, Fault_Theta_i, Fault_V_i, Fault_Friction_i, Fault_NormalStress, 
            Fault_V_Const, Fault_BulkIndex, FaultLengthStrike, FaultLengthDip, FaultStrikeAngle, 
            FaultDipAngle, FaultRakeAngle, FaultLengthStrike_Bulk, FaultLengthDip_Bulk, NormalStiffnessZero = 
                ReadSegmentInput(Input_Segment, FaultCount, RorT) 

            ############# Flip if needed and get Unit Vectors (triangle Only) ########
            if RorT == "T"
                P1, P2, P3, UnitVector_Normal, UnitVector_StrikeSlip, UnitVector_DipSlip, UnitVector_Slip = 
                    RotVerts_UnitVectors(Input_Segment, FaultCount, FaultRakeAngle) 
                    
            end


            if RorT == "R"
                ShearStiffness_H, NormalStiffness_H, Ranks_Shear, Ranks_Normal = 
                    HmatBuild_R(ShearModulus, PoissonRatio, ElementRange_SR,Input_Segment, Admissible, Tolerance)  
            elseif RorT == "T"
                ShearStiffness_H, NormalStiffness_H, Ranks_Shear, Ranks_Normal = 
                    HmatBuild_T(ShearModulus, PoissonRatio, ElementRange_SR, FaultRakeAngle, FaultCenter, UnitVector_Normal, 
                                UnitVector_Slip, Admissible, Tolerance, P1, P2, P3)   
            end


            save(OutputFileName, 
            "BlockI", BlockI,
            "BlockF", BlockF,
            "Ranks_Shear", Ranks_Shear,
            "Ranks_Normal", Ranks_Normal,
            "ShearStiffness_H",ShearStiffness_H, 
            "NormalStiffness_H", NormalStiffness_H)
            
            println("Saved File Name: ",OutputFileName)


        end
        Discritize()

        """

        fnamePart = "temp_Part$(DistributeIndex).jl"
        open(fnamePart, "w") do file
            write(file, CodeScript )
        end

    end



    if Sys.isapple()
        InputBulkFileNameCombine = "$(@__DIR__)/Input_BulkFaultGeometry.txt"
        OutputFileNameCombine = "$(@__DIR__)/Input_Discretized.jld2"
        HMatrixStructureFileCombine = "$(@__DIR__)/Input_HmatrixStructure.jld2"
    else
        InputBulkFileNameCombine ="Input_BulkFaultGeometry.txt"
        OutputFileNameCombine = "Input_Discretized.jld2"
        HMatrixStructureFileCombine = "Input_HmatrixStructure.jld2"
    end

    if Sys.isapple()
        LoadFileNameCombine = "$(@__DIR__)/temp_HMatPart_"
    else
        LoadFileNameCombine = "temp_HMatPart_"
    end

    ScriptCombine = """


        using DelimitedFiles
        using Base
        using PyPlot
        using PyCall
        using JLD2
        using LowRankApprox
        using Clustering
        using LinearAlgebra
        @pyimport matplotlib.patches as patches

        println("Merging parts. Do not close this window")

        pygui(true)

        include("scripts/Functions_BuildInputFile.jl")
        include("scripts/Functions_OKADA3D.jl")
        include("scripts/Functions_Hmatrix.jl")

        InputBulkFileName = "$(InputBulkFileNameCombine)"
        OutputFileName = "$(OutputFileNameCombine)"
        HMatrixStructureFile = "$(HMatrixStructureFileCombine)"

        Hmatrix = true
        ShearStiffness_H =[]
        NormalStiffness_H = []
        Ranks_Shear = []
        Ranks_Normal = []
        
        for Partindex=1:$(HowManyDistribution)

        LoadFileName = "$(LoadFileNameCombine)\$(Partindex).jld2"
        
        BlockI, BlockF, ShearStiffness_H_Part, NormalStiffness_H_Part, Ranks_Shear_Part, Ranks_Normal_Part = 
            load(LoadFileName, "BlockI", "BlockF", "ShearStiffness_H", "NormalStiffness_H",
            "Ranks_Shear", "Ranks_Normal")
            global ShearStiffness_H = [ShearStiffness_H ; ShearStiffness_H_Part]
            global NormalStiffness_H = [NormalStiffness_H ; NormalStiffness_H_Part]
            global Ranks_Shear = [Ranks_Shear ; Ranks_Shear_Part]
            global Ranks_Normal = [Ranks_Normal ; Ranks_Normal_Part]
        end


        ######################## Read General Input File #######################
        LoadingFaultCount, Admissible, ElementRange_SR, Switch_StrikeSlip_or_ReverseNormal, Input_Segment,
        ShearModulus, PoissonRatio, RockDensity,   MinimumNS, RorT, Input_Segment =
            load(HMatrixStructureFile, "LoadingFaultCount", "Admissible", "ElementRange_SR", "Switch_StrikeSlip_or_ReverseNormal", "Input_Segment",
                "ShearModulus", "PoissonRatio", "RockDensity",   "MinimumNS", "RorT", "Input_Segment")
        FaultCount = length(Input_Segment[:,1])

        ######################## Read Segmented Input File #######################
        FaultCenter, Fault_a, Fault_b, Fault_Dc, Fault_Theta_i, Fault_V_i, Fault_Friction_i, Fault_NormalStress, 
        Fault_V_Const, Fault_BulkIndex, FaultLengthStrike, FaultLengthDip, FaultStrikeAngle, 
        FaultDipAngle, FaultRakeAngle, FaultLengthStrike_Bulk, FaultLengthDip_Bulk, NormalStiffnessZero = 
            ReadSegmentInput(Input_Segment, FaultCount, RorT) 

        ############# Flip if needed and get Unit Vectors (triangle Only) ########
        if RorT == "T"
            P1, P2, P3, UnitVector_Normal, UnitVector_StrikeSlip, UnitVector_DipSlip, UnitVector_Slip = 
                RotVerts_UnitVectors(Input_Segment, FaultCount, FaultRakeAngle) 
        end



        save(OutputFileName, 
        "FaultCenter", FaultCenter,
        "ShearModulus", ShearModulus, "RockDensity", RockDensity, "PoissonRatio", PoissonRatio,
        "FaultLengthStrike", FaultLengthStrike, "FaultLengthDip", FaultLengthDip, "FaultStrikeAngle", FaultStrikeAngle, 
        "FaultDipAngle", FaultDipAngle, "FaultRakeAngle", FaultRakeAngle, "Fault_a", Fault_a, "Fault_b", Fault_b, "Fault_Dc", Fault_Dc, 
        "Fault_Theta_i", Fault_Theta_i, "Fault_V_i", Fault_V_i, "Fault_Friction_i", Fault_Friction_i, "Fault_NormalStress", Fault_NormalStress, 
        "Fault_V_Const", Fault_V_Const, "Fault_BulkIndex", Fault_BulkIndex, "FaultLengthStrike_Bulk", FaultLengthStrike_Bulk, 
        "FaultLengthDip_Bulk", FaultLengthDip_Bulk, "FaultCount", FaultCount, "LoadingFaultCount", LoadingFaultCount, 
        "Switch_StrikeSlip_or_ReverseNormal", Switch_StrikeSlip_or_ReverseNormal, "MinimumNormalStress", MinimumNS,
        "NormalStiffnessZero", NormalStiffnessZero, "RorT", RorT,   "Hmatrix", Hmatrix,
        "ShearStiffness_H", ShearStiffness_H, "NormalStiffness_H", NormalStiffness_H, "Admissible", Admissible, 
        "Ranks_Shear", Ranks_Shear, "Ranks_Normal", Ranks_Normal, "ElementRange_SR", ElementRange_SR    
        )
        println("Saved File Name: ",OutputFileName)

        ######################## Save Vertices for Triangles ##############################
        if RorT == "T"        
            file = jldopen(OutputFileName, "a+")
            write(file, "P1", P1) 
            write(file, "P2", P2) 
            write(file, "P3", P3) 
            close(file)
        end



    """

    fname = "temp_Combine.jl"
    open(fname, "w") do file
        write(file, ScriptCombine )
    end




    ######################## Paritial Discretization ##############################
    for DistributeIndex = 1 : HowManyDistribution    
        if isfile("temp_HMatPart_$DistributeIndex.jld2")
            rm("temp_HMatPart_$DistributeIndex.jld2")
        end
        # touch("DistributedDiscritization_beta/HMatPart_$DistributeIndex.jld2")
    end
    # touch("Input_Discretized.jld2")

    println("waiting for partial discretization")
    println("Progress can be monitored in separated windows")
    println("Do not close this REPL")


    for DistributeIndex = 1 : HowManyDistribution
        if Sys.islinux()
            run(`gnome-terminal -- bash -c "julia temp_Part$DistributeIndex.jl"`)
        elseif Sys.iswindows()
            run(`cmd /c start julia temp_Part$DistributeIndex.jl`)
        elseif Sys.isapple()

            program_to_run = "julia $(@__DIR__)/temp_Part$DistributeIndex.jl" 

            script = """
            tell application "Terminal"
                activate
                do script "$program_to_run"
            end tell
            """
            run(`osascript -e $script`)
        end
    end




    for DistributeIndex = 1 : HowManyDistribution
        while !isfile("temp_HMatPart_$DistributeIndex.jld2")
            sleep(1)
        end
    end

    ###############################################################################
    println("Combining the Parts")




    ######################## Combine Files ##############################
    if isfile("Input_Discretized.jld2")
        rm("Input_Discretized.jld2")
    end


    if Sys.islinux()
        run(`gnome-terminal -- bash -c "julia temp_Combine.jl"`)
    elseif Sys.iswindows()
        run(`cmd /c start julia temp_Combine.jl`)
    elseif Sys.isapple()           
        program_to_run = "julia $(@__DIR__)/temp_Combine.jl" 

        script = """
            tell application "Terminal"
            activate
            do script "$program_to_run"
            end tell
            """
        run(`osascript -e $script`)
    end


    for DistributeIndex = 1 : HowManyDistribution
        while !isfile("Input_Discretized.jld2")
            sleep(1)
        end
    end

    println("Discretization done")



    for DistributeIndex = 1 : HowManyDistribution
        rm("temp_Part$DistributeIndex.jl")
        rm("temp_HMatPart_$DistributeIndex.jld2")
    end
        rm("temp_Combine.jl")


end






if Hmatrix == true
    DistDisc()
elseif Hmatrix == false
    DiscritizeFullMatrix()
end




