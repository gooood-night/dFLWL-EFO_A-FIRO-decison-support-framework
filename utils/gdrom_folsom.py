# this file includes the pre-trained GDROM operation rules for Folsom Lake
# For Folsom Lake, GDROM identified two representative operation modules: one module (M1) mainly for normal operations and another module (M0) mainly for flood control, along with their detailed transition conditions. 

#******************************
def GDROM_M1(I1, S0):
    Inflow = I1 * 1000 # TAF to AF
    Storage = S0 * 1000
    
    if (Storage <= 673955.0):
        Release = 3932.8
    if (Storage > 673955.0):
        Release = 6777.0
    
    R1 = Release/1000
    return R1


#******************************
def GDROM_M0(I1, S0):
    Inflow = I1 * 1000 # TAF to AF
    Storage = S0 * 1000
    
    if (Inflow <= 28060.7) and (Inflow <= 13723.0) and (Storage <= 791778.5): 
        Release = 9081.2
    if (Inflow <= 28060.7) and (Inflow <= 13723.0) and (Storage > 791778.5): 
        Release = 13525.6
    if (Inflow <= 28060.7) and (Inflow > 13723.0) and (Inflow <= 20384.8) and (Storage <= 564218.5): 
        Release = 20066.6
    if (Inflow <= 28060.7) and (Inflow > 13723.0) and (Inflow <= 20384.8) and (Storage > 564218.5): 
        Release = 14621.6
    if (Inflow <= 28060.7) and (Inflow > 13723.0) and (Inflow > 20384.8) and (Storage <= 636247.5): 
        Release = 26320.1
    if (Inflow <= 28060.7) and (Inflow > 13723.0) and (Inflow > 20384.8) and (Storage > 636247.5): 
        Release = 19352.3
    if (Inflow > 28060.7) and (Inflow <= 62442.0) and (Inflow <= 38817.4): 
        Release = 31473.6
    if (Inflow > 28060.7) and (Inflow <= 62442.0) and (Inflow > 38817.4): 
        Release = 43771.8
    if (Inflow > 28060.7) and (Inflow > 62442.0) and (Storage <= 660068.0): 
        Release = 56287.8
    if (Inflow > 28060.7) and (Inflow > 62442.0) and (Storage > 660068.0): 
        Release = 102071.1
    
    R1 = Release/1000
    return R1


#******************************
def GDROM_module_transition(I1, S0, PDSI):
    Inflow = I1 * 1000 # TAF to AF
    Storage = S0 * 1000 # TAF to TA
    
    if (Inflow <= 9661.1): 
        module = 1
    if (Inflow > 9661.1) and (Storage <= 640763.0) and (PDSI <= -2.18): 
        module = 1
    if (Inflow > 9661.1) and (Storage > 640763.0) and (PDSI <= -2.18) and (DOY <= 82.0): 
        module = 0
    if (Inflow > 9661.1) and (Storage > 640763.0) and (PDSI <= -2.18) and (DOY > 82.0): 
        module = 1
    if (Inflow > 9661.1) and (PDSI <= 2.3) and (PDSI > -2.18): 
        module = 0
    if (Inflow > 9661.1) and (PDSI <= 2.4) and (PDSI > 2.3): 
        module = 1
    if (Inflow > 9661.1) and (Storage <= 346370.0) and (PDSI > 2.4): 
        module = 1
    if (Inflow > 9661.1) and (Storage > 346370.0) and (PDSI > 2.4): 
        module = 0
        
    return module


#******************************
def GDROM_release_prediction(I1, S0, PDSI_obs):
    
    #%%%%%%%  using the GDROM for non-flood season release decisions
    print("GDROM")
    module_num = GDROM_module_transition(I1, S0, PDSI_obs)
    if module_num == 0:
        R1 =  GDROM_M0(I1, S0)
    else:
        R1 =  GDROM_M1(I1, S0)
        
    return R1