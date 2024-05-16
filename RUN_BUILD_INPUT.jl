

using DelimitedFiles
using Base
using PyPlot
using PyCall
using JLD2
pygui(true)

include("Functions_BuildInputFile.jl")
include("Functions_OKADA3D.jl")

function BuildInputFromBulkGeometry()


    InputBulkFileName="Input_BulkFaultGeometry.txt"
    OutputFileName="Input_Discretized.jld2"


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
    for i=1:size(Input_Bulk,1)
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



    ########################## Build Stiffness Matrix ##############################

    ElementPartRoughCount = 2000
    StiffnessMatrixShearOriginal=zeros(FaultCount,FaultCount)
    StiffnessMatrixNormalOriginal=zeros(FaultCount,FaultCount)

    if Switch_StrikeSlip_or_ReverseNormal == 1
        println("Preparing for discretization")
        StiffnessMatrixShearOriginal, StiffnessMatrixNormalOriginal = 
            BuildMatrixByPartsShear(FaultCount, ElementPartRoughCount, Input_Segment,  ShearModulus, PoissonRatio)
        
    elseif Switch_StrikeSlip_or_ReverseNormal == 2
        println("Preparing for discretization")
        StiffnessMatrixShearOriginal, StiffnessMatrixNormalOriginal = 
            BuildMatrixByPartsNormal(FaultCount, ElementPartRoughCount, Input_Segment,  ShearModulus, PoissonRatio)

    else
        println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        println("!!!!!!! Switch_StrikeSlip_or_ReverseNormal should be either 1 or 2  !!!!!!!")
        println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    end



    println("Stiffness NaN count: ",sum(isnan.(StiffnessMatrixShearOriginal)))

    ########################## Remove Unstable Faults ##############################    
    ReducedStiffnessMatrixShear, ReducedStiffnessMatrixNormal, ReducedInput_Segment=
    CheckTooClose(StiffnessMatrixShearOriginal, StiffnessMatrixNormalOriginal, Input_Segment, Input_Bulk, DropCrit, DropCritNormalStressMultiplier);


    ################################### Save Files #################################
    SaveResults(ReducedStiffnessMatrixShear, ReducedStiffnessMatrixNormal, ReducedInput_Segment,
         OutputFileName, ShearModulus, PoissonRatio, RockDensity, Switch_StrikeSlip_or_ReverseNormal, MinimumNS);

end




Input = BuildInputFromBulkGeometry()