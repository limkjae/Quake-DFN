###################################################################################
######################## input parameters for the Simulation ######################


######################### Simulation Time Set #############################
TotalStep = 5000 # Total simulation step
SaveStep = 5000 # Automatically saved every this step
RecordStep = 10 # Simulation sampling rate !! should be a factor of SaveStep !!

ThreadCount = 5 # Only used for HMatrix Simulations. 5~8 seems to be most efficient in most of the computer

########################## Time Stepping Setup ############################
DtCut = 5
SwitchV = 1e-2
RuptureTimeStepMultiple = 3
MaximumDt = 1e7
VerticalLengthScaleforM = 0 # if 0, Mass is automatically determined based on the fault length (radiation damping dominated for large rupture). If not, M = VerticalLengthScaleforM * density / 2


##########################  Save File Name  ###############################
SaveFileName = "Result"


####### Strong Interaction Supression for Numerical Stability in case too Strong Interaction ######
### If interaction is larger than "StrongInteractionCriteriaMultiple" * self stiffness,
### the interaction will be reduced to "StrongInteractionCriteriaMultiple" * self stiffness
### Only applied when larger than 0. The higher, the more tolerance of strong interaction. 
StrongInteractionCriteriaMultiple = 0







if Sys.islinux()
    run(`gnome-terminal -- bash -c "julia --threads $(ThreadCount) scripts/QUAKEDFN_RunARG.jl $(TotalStep) $(SaveStep) $(RecordStep) $(ThreadCount) $(DtCut) $(SwitchV) $(RuptureTimeStepMultiple) $(MaximumDt) $(VerticalLengthScaleforM) $(SaveFileName) $(StrongInteractionCriteriaMultiple)"`)
elseif Sys.iswindows()
    run(`cmd /c start julia --threads $(ThreadCount) scripts/QUAKEDFN_RunARG.jl $(TotalStep) $(SaveStep) $(RecordStep) $(ThreadCount) $(DtCut) $(SwitchV) $(RuptureTimeStepMultiple) $(MaximumDt) $(VerticalLengthScaleforM) $(SaveFileName) $(StrongInteractionCriteriaMultiple)`)

elseif Sys.isapple()
        program_to_run = "julia --threads $(ThreadCount) $(@__DIR__)/scripts/QUAKEDFN_RunARG.jl $(TotalStep) $(SaveStep) $(RecordStep) $(ThreadCount) $(DtCut) $(SwitchV) $(RuptureTimeStepMultiple) $(MaximumDt) $(VerticalLengthScaleforM) $(SaveFileName) $(StrongInteractionCriteriaMultiple)"
        script = """
            tell application "Terminal"
            activate
            do script "$program_to_run"
            end tell
            """
        run(`osascript -e $script`)
    # run(`open -a Terminal julia --threads $(ThreadCount) scripts/QUAKEDFN_RunARG.jl $(TotalStep) $(SaveStep) $(RecordStep) $(ThreadCount) $(DtCut) $(SwitchV) $(RuptureTimeStepMultiple) $(MaximumDt) $(VerticalLengthScaleforM) $(SaveFileName) $(StrongInteractionCriteriaMultiple)`)
end

println("Simulation starts in a new terminal")