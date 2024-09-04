

using DelimitedFiles
using Base
using PyPlot
using PyCall
using JLD2
using LowRankApprox
using Clustering
using LinearAlgebra
using Distributed
@pyimport matplotlib.patches as patches

pygui(true)

include("Functions_BuildInputFile.jl")
include("Functions_OKADA3D.jl")
include("Results/Functions_Plot.jl")
include("Functions_Hmatrix.jl")


function BuildInputFromBulkGeometry_H()

    InputBulkFileName="Input_BulkFaultGeometry.txt"
    OutputFileName="Input_Discretized.jld2"
    HMatrixStructureFile = "Input_HmatrixStructure.jld2"
    ##########################################################################
    ########################## Hmatrix compress? #############################
    HMatrixCompress = 1 # If this is 1, stiffness Matrix will be compressed using Input_HmatrixStructure.jld2
    SaveOriginalMatrix = 0 # 1: save Original Matrix (can be very large), 0: Discard Original Matrix. 
    
    #####----- Hmatrix compression Tolerance ----#####
    Tolerance = 1e3 # pascal for 1m slip (More approximaion for higher Tolerance). 

    ##########################################################################
    ##########################################################################

    #####---------   Plots?  --------#####
    Plot_HMatStructure = 0 # HMatrix structure plot
    Plot3D_HMatrixBlock = 0 # Blocks in 3D plot
    Plot3D_DiscretizedBlock = 0 # 3D discritized block
    
    ######################################

    ElementPartRoughCount = 2000
    

    if HMatrixCompress == 1
        if isfile(HMatrixStructureFile)            
            println("H Matrix will be built")
            Block_Ctr_Diam, Block_Range_Level, Input_Segment, LoadingFaultExist, HowManyDivisionEachLevel,
            MinimumElementsToCut, ArrangePoint, Admissible, ElementRange_SR, Switch_StrikeSlip_or_ReverseNormal, 
            ShearModulus, PoissonRatio, RockDensity, DropCrit, DropCritNormalStressMultiplier, MinimumNS, Input_Bulk =
                load(HMatrixStructureFile, "Block_Ctr_Diam", "Block_Range_Level", "Input_Segment", "LoadingFaultExist", "HowManyDivisionEachLevel",
                     "MinimumElementsToCut", "ArrangePoint", "Admissible", "ElementRange_SR", "Switch_StrikeSlip_or_ReverseNormal", 
                     "ShearModulus", "PoissonRatio", "RockDensity", "DropCrit", "DropCritNormalStressMultiplier", "MinimumNS", "Input_Bulk")

            Input_BulkfromFile=readdlm(InputBulkFileName)
            if Input_Bulk != Input_BulkfromFile
                println("???????????????????????????????????????????????????????????????????????????????????????????")
                println("???          H Matrix Structure File dosen't correspond to Bulk Input File              ???")
                println("???????????????????????????????????????????????????????????????????????????????????????????")
            end
            FaultCount = length(Input_Segment[:,1])


        else 
            # println("H Matrix Structure File does not exist")
            println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            println("!!!                        H Matrix Structure File does not exist                       !!!")
            println("!!!      To use Hmatrix apploximation, Run 'RUN_BUILD_HMatrix_Structure.jl' first       !!!")
            println("!!!                          Otherwise, set  HMatrixCompress = 0                        !!!")
            println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            error()

        end
    else
        ############################# Read Bulk Input ##################################
        ######++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++######
        ## Bulk File Order 
        ##  1.SwitchSSRN  2.ShearMod  3.PoiRat  4.R_Density   
        ##  5. Crit_TooClose   6. TooCloseNormal_Multiplier
        ##  ----------------------------------------------------------------------------
        ##  1.Ctr_X     2.Ctr_Y 3.Ctr_Z 4.St_L	    5.Dip_L	    6.StAng	    7.DipAng	8.LR/RN
        ##  9.a         10.b	11.Dc	12.Theta_i	13. V_i     14. Friction_i 15.NormalStress at surface [Pa]  
        ##  16. NoarmalStress Gradient [Pa] 17. V_Const     18. Minimum Segment Length

        Input_Bulk=readdlm(InputBulkFileName)
        
        Switch_StrikeSlip_or_ReverseNormal = Input_Bulk[2,1] 
        ShearModulus = Input_Bulk[2,2]
        PoissonRatio = Input_Bulk[2,3]
        RockDensity = Input_Bulk[2,4]
        DropCrit= Input_Bulk[2,5]
        DropCritNormalStressMultiplier= Input_Bulk[2,6]
        MinimumNS=Input_Bulk[2,7]

        Input_Bulk=Input_Bulk[4:end,:]
        Input_Bulk=Input_Bulk[sortperm(Input_Bulk[:, 17]), :]
        

        # Adjust if positive depth exists
        for i in eachindex(Input_Bulk[:,1])
            if Input_Bulk[i,3]<Input_Bulk[i,5]/2*sind(Input_Bulk[i,7])
                println("Caution! Fault ",i," may have negative depth")
            end
        end

        ########^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^########
        ################################################################################


        ############################### Bulk to Segment ################################
        ######++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++######
        ## Segment File Order 
        ##  1.Ctr_X     2.Ctr_Y 3.Ctr_Z 4.St_L	    5.Dip_L	    6.StAng	    7.DipAng	8.LR/RN
        ##  9.a         10.b	11.Dc	12.Theta_i	13. V_i     14. Friction_i 15.NormalStress  
        ##  16. V_Const 17. Bulk Index     18. Bulk Strike Length      19. Bulk Dip Length
        Input_Segment = BulkToSegment(Input_Bulk);
        FaultCount=   size(Input_Segment,1)
        Input_Segment = [Input_Segment ones(FaultCount) zeros(FaultCount)]
        Input_Segment = Input_Segment[sortperm(Input_Segment[:, 16]), :] # move the loading faults to the end
    end


    if HMatrixCompress == 1
        #####---------------------------    Hmat Plots  -------------------------------#####
        #####                               3D Plot                                    #####
        if Plot3D_HMatrixBlock == 1
            println("Plotting 3D Group")
            PlotRotation=[45,-30]
            Edge = 0
            Transparent = 0
            MinMax_Axis=0 # automatically detect max and min 
            LoadingFaultCount=0 
            ColorMinMax = 0  
            figure(1)
            clf()
            MaxVaule, MinValue = FaultPlot_3D_Color_General(Input_Segment[:,1:3],Input_Segment[:,4], Input_Segment[:,5],
                Input_Segment[:,6], Input_Segment[:,7], Input_Segment[:,8], Input_Segment[:,20], 
                PlotRotation, MinMax_Axis, ColorMinMax, Transparent, Edge, LoadingFaultCount)

            plotforcbar=  scatter([1,1],[1,1],0.1, [MinValue,MaxVaule], cmap="jet")
            colorbar(plotforcbar, pad=0.15)
            figure(1).canvas.draw()
            xlabel("x")
            ylabel("y")
        end

        #####                     Hierarchical Matrix Structure Plot                  #####
        if Plot_HMatStructure == 1
            figure(2)
            println("Plotting Hmatrix Structure")
            clf()
            ax = gca()
            # ax[:set_aspect]("equal")
            for i=1:length(ElementRange_SR[:,1])
                if Admissible[i] > 0   
                    c = PyObject(patches.Rectangle((ElementRange_SR[i,3]-1, -ElementRange_SR[i,1]+1), ElementRange_SR[i,4] - ElementRange_SR[i,3]+1, 
                                -ElementRange_SR[i,2] + ElementRange_SR[i,1]-1, linewidth=1, edgecolor="k", facecolor=[0.4  0.4  1]))    
                else           
                    c = PyObject(patches.Rectangle((ElementRange_SR[i,3]-1, -ElementRange_SR[i,1]+1), ElementRange_SR[i,4] - ElementRange_SR[i,3]+1, 
                                -ElementRange_SR[i,2] + ElementRange_SR[i,1]-1, linewidth=1, edgecolor="k", facecolor=[1 0.4 0.4]))
                end
                # ax.text( (ElementRange_SR[i,3] + ElementRange_SR[i,4])/2, -(ElementRange_SR[i,1] + ElementRange_SR[i,2])/2, i ,size=8, horizontalalignment="center", verticalalignment="center", color="k") 
                ax.add_patch(c) 
            end
            xlim(0,FaultCount)
            ylim(-FaultCount,0)
        end
        ########^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^########
    end
    ############################ end of Hmatrix Sorting ################################



    ####################################################################################
    ############################# Build Stiffness Matrix ###############################
    if HMatrixCompress == 1
        if SaveOriginalMatrix == 0
            BlockCount = length(ElementRange_SR[:,1])
            ShearStiffness_H = Any[0]
            NormalStiffness_H = Any[0]
            Ranks_Shear = zeros(Int, BlockCount)
            Ranks_Normal = zeros(Int, BlockCount)

            TotalElments = FaultCount * FaultCount 
            println("Building Hmatrix Block by Block")
            println("Full Matrix Will not be saved")
            println("Preparing for discretization")
            BlockSize = 0.0
            
            # @sync Threads.@threads for BlockIndex = 1: BlockCount
            for BlockIndex = 1: BlockCount
                
                Input_SegmentS = Input_Segment[ElementRange_SR[BlockIndex,1]:ElementRange_SR[BlockIndex,2],:]
                Input_SegmentR = Input_Segment[ElementRange_SR[BlockIndex,3]:ElementRange_SR[BlockIndex,4],:]
                SourceCount = length(Input_SegmentS[:,1])
                ReceiverCount = length(Input_SegmentR[:,1])
                StiffnessMatrixShearThisBlock=zeros(ReceiverCount,SourceCount)
                StiffnessMatrixNormalThisBlock=zeros(ReceiverCount,SourceCount)

                DivisionCountS = round(Int,SourceCount / ElementPartRoughCount)
                DivisionCountR = round(Int,ReceiverCount / ElementPartRoughCount)
                if DivisionCountS == 0; DivisionCountS =1; end
                if DivisionCountR == 0; DivisionCountR =1; end
                PartedElementCountS = SourceCount ÷ DivisionCountS
                PartedElementCountR = ReceiverCount ÷ DivisionCountR
                TotalParts = DivisionCountS * DivisionCountR
                CurrentPart = 0
                for i=1:DivisionCountS
                    for j=1:DivisionCountR
                        CurrentPart =  CurrentPart +1
                        Init_S = (i-1)*PartedElementCountS + 1
                        Fin_S = i*PartedElementCountS
                        Init_R =  (j-1)*PartedElementCountR + 1
                        Fin_R = j*PartedElementCountR
                        if i == DivisionCountS; Fin_S = SourceCount; end
                        if j == DivisionCountR; Fin_R = ReceiverCount; end
                        BlockSize = BlockSize + (Fin_S - Init_S) * (Fin_R - Init_R)
                        println("Part: ", CurrentPart,"/",DivisionCountS*DivisionCountR, " BlockIndex: ",BlockIndex, "/",BlockCount," Progress:",BlockSize/TotalElments)

                        if Switch_StrikeSlip_or_ReverseNormal == 1
                        StiffnessMatrixShearThisBlock[Init_R:Fin_R,Init_S:Fin_S], StiffnessMatrixNormalThisBlock[Init_R:Fin_R,Init_S:Fin_S] = 
                        StiffnessMatrix_ByParts_Calculation_StrikeSlip(Input_SegmentS[Init_S:Fin_S,:], Input_SegmentR[Init_R:Fin_R,:], ShearModulus, PoissonRatio,
                                                            CurrentPart, TotalParts)

                        elseif Switch_StrikeSlip_or_ReverseNormal == 2                                                        
                        StiffnessMatrixShearThisBlock[Init_R:Fin_R,Init_S:Fin_S], StiffnessMatrixNormalThisBlock[Init_R:Fin_R,Init_S:Fin_S] = 
                        StiffnessMatrix_ByParts_Calculation_NormalReverse(Input_SegmentS[Init_S:Fin_S,:], Input_SegmentR[Init_R:Fin_R,:], ShearModulus, PoissonRatio,
                                                            CurrentPart, TotalParts)                                                                    
                        else
                            println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                            println("!!!!!!! Switch_StrikeSlip_or_ReverseNormal should be either 1 or 2  !!!!!!!")
                            println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                        end
                    end
                end                        
                
                ######################################################################
                ############################# Compress ###############################
                if Admissible[BlockIndex] > 0
                    ApproxMatrixS = pqrfact(StiffnessMatrixShearThisBlock, atol = Tolerance)
                    push!(ShearStiffness_H,ApproxMatrixS)
                    Ranks_Shear[BlockIndex] = size(ApproxMatrixS[:Q],2)
                    
                    ApproxMatrixN = pqrfact(StiffnessMatrixNormalThisBlock, atol = Tolerance)
                    push!(NormalStiffness_H,ApproxMatrixN)
                    Ranks_Normal[BlockIndex] = size(ApproxMatrixN[:Q],2)
                else 
                    push!(ShearStiffness_H,StiffnessMatrixShearThisBlock)
                    push!(NormalStiffness_H,StiffnessMatrixNormalThisBlock)
                end
            end
            ShearStiffness_H = ShearStiffness_H[2:end]
            NormalStiffness_H = NormalStiffness_H[2:end]
            
        else 
            ################### Original Matrix Should be Saved #####################
            println("Build Full Matrix and then compress")

            StiffnessMatrixShearOriginal=zeros(FaultCount,FaultCount)
            StiffnessMatrixNormalOriginal=zeros(FaultCount,FaultCount)    
            if Switch_StrikeSlip_or_ReverseNormal == 1

                println("Preparing for discretization")
                StiffnessMatrixShearOriginal, StiffnessMatrixNormalOriginal = 
                    BuildMatrixByPartsStrikeSlip(FaultCount, ElementPartRoughCount, Input_Segment,  ShearModulus, PoissonRatio)
                
            elseif Switch_StrikeSlip_or_ReverseNormal == 2
                println("Preparing for discretization")
                StiffnessMatrixShearOriginal, StiffnessMatrixNormalOriginal = 
                    BuildMatrixByPartsNormalReverse(FaultCount, ElementPartRoughCount, Input_Segment,  ShearModulus, PoissonRatio)
    
            else
                println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                println("!!!!!!! Switch_StrikeSlip_or_ReverseNormal should be either 1 or 2  !!!!!!!")
                println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            end

            BlockCount = length(ElementRange_SR[:,1])
            ShearStiffness_H = Any[0]
            Ranks_Shear = zeros(Int, BlockCount)
            NormalStiffness_H = Any[0]
            Ranks_Normal = zeros(Int, BlockCount)
            BlockIndex = 0
            println("compressing")
            for i=1:BlockCount
                BlockIndex = BlockIndex + 1
                
                if Admissible[BlockIndex] > 0
                    OrigianlMatrixToApproximate = StiffnessMatrixShearOriginal[ElementRange_SR[i,3]:ElementRange_SR[i,4],ElementRange_SR[i,1]:ElementRange_SR[i,2]]
                    ApproxMatrixS = pqrfact(OrigianlMatrixToApproximate, atol = Tolerance)
                    push!(ShearStiffness_H,ApproxMatrixS)
                    Ranks_Shear[BlockIndex] = size(ApproxMatrixS[:Q],2)
                    
                    OrigianlMatrixToApproximate = StiffnessMatrixNormalOriginal[ElementRange_SR[i,3]:ElementRange_SR[i,4],ElementRange_SR[i,1]:ElementRange_SR[i,2]]
                    ApproxMatrixN = pqrfact(OrigianlMatrixToApproximate, atol = Tolerance)
                    push!(NormalStiffness_H,ApproxMatrixN)
                    Ranks_Normal[BlockIndex] = size(ApproxMatrixN[:Q],2)
                else 
                    push!(ShearStiffness_H,StiffnessMatrixShearOriginal[ElementRange_SR[i,3]:ElementRange_SR[i,4],ElementRange_SR[i,1]:ElementRange_SR[i,2]])
                    push!(NormalStiffness_H,StiffnessMatrixNormalOriginal[ElementRange_SR[i,3]:ElementRange_SR[i,4],ElementRange_SR[i,1]:ElementRange_SR[i,2]])
                end
    
            end
            ShearStiffness_H = ShearStiffness_H[2:end]
            NormalStiffness_H = NormalStiffness_H[2:end]

        end


    else
        ################## No Hmatrix Compression ##################
        println("building full stiffness matrix (No-Hmatrix Compression)")
        StiffnessMatrixShearOriginal=zeros(FaultCount,FaultCount)
        StiffnessMatrixNormalOriginal=zeros(FaultCount,FaultCount)

        if Switch_StrikeSlip_or_ReverseNormal == 1
            println("Preparing for discretization")
            StiffnessMatrixShearOriginal, StiffnessMatrixNormalOriginal = 
                BuildMatrixByPartsStrikeSlip(FaultCount, ElementPartRoughCount, Input_Segment,  ShearModulus, PoissonRatio)
            
        elseif Switch_StrikeSlip_or_ReverseNormal == 2
            println("Preparing for discretization")
            StiffnessMatrixShearOriginal, StiffnessMatrixNormalOriginal = 
                BuildMatrixByPartsNormalReverse(FaultCount, ElementPartRoughCount, Input_Segment,  ShearModulus, PoissonRatio)

        else
            println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            println("!!!!!!! Switch_StrikeSlip_or_ReverseNormal should be either 1 or 2  !!!!!!!")
            println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        end

    end



    NormalStiffnessZero = 0    

    ################################################################################
    ################################### Save Files #################################

    FaultCenter = Input_Segment[:,1:3]
    FaultLengthStrike = Input_Segment[:,4]
    FaultLengthDip = Input_Segment[:,5]
    FaultStrikeAngle = Input_Segment[:,6]
    FaultDipAngle = Input_Segment[:,7]
    FaultLLRR = Input_Segment[:,8]
    Fault_a = Input_Segment[:,9]
    Fault_b = Input_Segment[:,10]
    Fault_Dc = Input_Segment[:,11]
    Fault_Theta_i = Input_Segment[:,12]
    Fault_V_i = Input_Segment[:,13]
    Fault_Friction_i = Input_Segment[:,14]
    Fault_NormalStress = Input_Segment[:,15]
    Fault_V_Const = Input_Segment[:,16]
    Fault_BulkIndex = Input_Segment[:,17]
    FaultLengthStrike_Bulk = Input_Segment[:,18]
    FaultLengthDip_Bulk = Input_Segment[:,19]

    LoadingFaultCount=length(Fault_V_Const[Fault_V_Const.>0])

    if HMatrixCompress ==1
        if SaveOriginalMatrix ==1 
            save(OutputFileName, 
            "StiffnessMatrixShear", StiffnessMatrixShearOriginal, "StiffnessMatrixNormal", StiffnessMatrixNormalOriginal, "FaultCenter", FaultCenter,
            "ShearModulus", ShearModulus, "RockDensity", RockDensity, "PoissonRatio", PoissonRatio,
            "FaultLengthStrike", FaultLengthStrike, "FaultLengthDip", FaultLengthDip, "FaultStrikeAngle", FaultStrikeAngle, 
            "FaultDipAngle", FaultDipAngle, "FaultLLRR", FaultLLRR, "Fault_a", Fault_a, "Fault_b", Fault_b, "Fault_Dc", Fault_Dc, 
            "Fault_Theta_i", Fault_Theta_i, "Fault_V_i", Fault_V_i, "Fault_Friction_i", Fault_Friction_i, "Fault_NormalStress", Fault_NormalStress, 
            "Fault_V_Const", Fault_V_Const, "Fault_BulkIndex", Fault_BulkIndex, "FaultLengthStrike_Bulk", FaultLengthStrike_Bulk, 
            "FaultLengthDip_Bulk", FaultLengthDip_Bulk, "FaultCount", FaultCount, "LoadingFaultCount", LoadingFaultCount, 
            "Switch_StrikeSlip_or_ReverseNormal", Switch_StrikeSlip_or_ReverseNormal, "MinimumNormalStress", MinimumNS,
            "Ranks_Shear", Ranks_Shear, "Ranks_Normal",Ranks_Normal,"ElementRange_SR", ElementRange_SR, "ShearStiffness_H",ShearStiffness_H, "NormalStiffness_H", NormalStiffness_H, "Admissible", Admissible,
            "NormalStiffnessZero", NormalStiffnessZero,"SaveOriginalMatrix",SaveOriginalMatrix)
            println("Saved File Name: ",OutputFileName)
        else
            save(OutputFileName, 
            "FaultCenter", FaultCenter,
            "ShearModulus", ShearModulus, "RockDensity", RockDensity, "PoissonRatio", PoissonRatio,
            "FaultLengthStrike", FaultLengthStrike, "FaultLengthDip", FaultLengthDip, "FaultStrikeAngle", FaultStrikeAngle, 
            "FaultDipAngle", FaultDipAngle, "FaultLLRR", FaultLLRR, "Fault_a", Fault_a, "Fault_b", Fault_b, "Fault_Dc", Fault_Dc, 
            "Fault_Theta_i", Fault_Theta_i, "Fault_V_i", Fault_V_i, "Fault_Friction_i", Fault_Friction_i, "Fault_NormalStress", Fault_NormalStress, 
            "Fault_V_Const", Fault_V_Const, "Fault_BulkIndex", Fault_BulkIndex, "FaultLengthStrike_Bulk", FaultLengthStrike_Bulk, 
            "FaultLengthDip_Bulk", FaultLengthDip_Bulk, "FaultCount", FaultCount, "LoadingFaultCount", LoadingFaultCount, 
            "Switch_StrikeSlip_or_ReverseNormal", Switch_StrikeSlip_or_ReverseNormal, "MinimumNormalStress", MinimumNS,
            "Ranks_Shear", Ranks_Shear, "Ranks_Normal",Ranks_Normal,"ElementRange_SR", ElementRange_SR, "ShearStiffness_H",ShearStiffness_H, "NormalStiffness_H", NormalStiffness_H, "Admissible", Admissible,
            "NormalStiffnessZero", NormalStiffnessZero,"SaveOriginalMatrix",SaveOriginalMatrix)
            println("Saved File Name: ",OutputFileName)
        end

    else 
            save(OutputFileName, 
            "StiffnessMatrixShear", StiffnessMatrixShearOriginal, "StiffnessMatrixNormal", StiffnessMatrixNormalOriginal, "FaultCenter", FaultCenter,
            "ShearModulus", ShearModulus, "RockDensity", RockDensity, "PoissonRatio", PoissonRatio,
            "FaultLengthStrike", FaultLengthStrike, "FaultLengthDip", FaultLengthDip, "FaultStrikeAngle", FaultStrikeAngle, 
            "FaultDipAngle", FaultDipAngle, "FaultLLRR", FaultLLRR, "Fault_a", Fault_a, "Fault_b", Fault_b, "Fault_Dc", Fault_Dc, 
            "Fault_Theta_i", Fault_Theta_i, "Fault_V_i", Fault_V_i, "Fault_Friction_i", Fault_Friction_i, "Fault_NormalStress", Fault_NormalStress, 
            "Fault_V_Const", Fault_V_Const, "Fault_BulkIndex", Fault_BulkIndex, "FaultLengthStrike_Bulk", FaultLengthStrike_Bulk, 
            "FaultLengthDip_Bulk", FaultLengthDip_Bulk, "FaultCount", FaultCount, "LoadingFaultCount", LoadingFaultCount,
            "Switch_StrikeSlip_or_ReverseNormal", Switch_StrikeSlip_or_ReverseNormal, "MinimumNormalStress", MinimumNS,
            "NormalStiffnessZero", NormalStiffnessZero)
            println("Saved File Name: ",OutputFileName)
    
    end


    if Plot3D_DiscretizedBlock == 1
        figure(3)
        clf()
        FaultPlot_3D(Input_Segment[:,1:3],Input_Segment[:,4], Input_Segment[:,5], 
                    Input_Segment[:,6], Input_Segment[:,7], Input_Segment[:,8])
            xlabel("x")
            ylabel("y")
        figure(3).canvas.draw()
    end




end


Input = BuildInputFromBulkGeometry_H()
