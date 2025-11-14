###################################################################################
######################## input parameters for the Simulation ######################


######################### Simulation Time Set #############################
TotalStep = 2000 # Total simulation step
SaveStep = 1000 # Automatically saved every this step
RecordStep = 10 # Simulation sampling rate !! should be a factor of SaveStep !!

ThreadCount = 5 # Only used in HMatrix Simulations

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







#################### Run simulation in a new terminal #####################
 run(`cmd /c start julia --threads $(ThreadCount) scripts/QUAKEDFN_RunARG.jl 
    $(TotalStep) $(SaveStep) $(RecordStep) $(ThreadCount) $(DtCut) 
    $(SwitchV) $(RuptureTimeStepMultiple) $(MaximumDt) $(VerticalLengthScaleforM) $(SaveFileName) $(StrongInteractionCriteriaMultiple)`)

println("Simulation starts in a new terminal")