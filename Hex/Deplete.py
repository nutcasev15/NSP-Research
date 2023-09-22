# OpenMC 0.13.3 Python Master Script for MF Ratio Depletion Study


############### Library Imports
import os
import json
import pickle
import numpy as np
import openmc
import openmc.deplete


############### Import Model Builder Function
from OMC_model import OMC_model


############### Define Investigation Parameters
# Define Save File Name
Pikname = 'Dep-Ref'

# Define Maximum Fuel Element Linear Power
pow = 275 # W / cm

# Define List of 8 MF Volume Ratios between 1 and 25
MF_vol = np.logspace(0, np.log10(25), 8).tolist()

# Define Depletion Schedule up to 10 Years
dep_sch = np.logspace(-3, -1, 3) # Setup Schedule to Equilibriate Fission Products 
dep_sch = np.append(dep_sch, np.logspace(0, 1, 3)).tolist() # Deplete Until 10 a
dep_steps = np.diff(np.array(dep_sch), prepend=0).tolist() # Get Timesteps for CELI Integrator

# Define Depletion Chain File Name for OpenMC
chain = 'chain_casl_pwr'

# Define List to Store Moderator Name from Model Builder Function
Mname = []

# Define List to Store MF Density Ratio from Model Builder Function
MF_den_ratio = []

# Define List to Store Coolant Channel Diameter from Model Builder Function
chnl_dia = []

# Define List to Store Fuel Element Circumscribing Diameter from Model Builder Function
hex_dia = []

# Define List to Store Fuel Element Linear Mass from Model Builder Function
hex_mass = []

# Define Lists to Store Kinf Results
K_EVL = []
K_BOL = []
K_EOL = []

# Define List to Store U Atom Burnup
U_burn = []


############### Import Serpent Fission Q Data for Depletion Simulations
with open('serpent_fissq.json', 'r') as f:
    serpent_fission_q = json.load(f)


############### Run Depletion Routine for Each Moderator
for ID in range(0, 7):
    # Define Temporary Containers for Model Builder and Calculation Outputs
    Modname = []
    MF_dens = []
    cc_dia = []
    hp_dia = []
    fe_mass = []
    kinf = []
    burn = []

    for MF in MF_vol:
       ######## Call Model Builder Function
        Model_Data = OMC_model(ID, MF)

        ######## Save Model Output Data
        Modname.append(Model_Data[1])
        MF_dens.append(Model_Data[2])
        cc_dia.append(Model_Data[3])
        hp_dia.append(Model_Data[4])
        fe_mass.append(Model_Data[5])

        ######## Clear Run Directory of Previous Results
        os.system('rm -rf *.xml')
        os.system('rm -rf *.h5')

        ######## Override Simulation Settings to Reduce Runtime
        Model_Data[0].settings.batches = 15
        Model_Data[0].settings.inactive = 5
        Model_Data[0].settings.particles = 5000

        ######## Initialise OpenMC's Coupled Depletion Operator with a CELI Integrator
        dep_opr = openmc.deplete.CoupledOperator(Model_Data[0],
                                                 chain,
                                                 fission_yield_mode='average',
                                                 fission_q=serpent_fission_q,
                                                 diff_burnable_mats=True)
        dep_int = openmc.deplete.CELIIntegrator(dep_opr,
                                                dep_steps,
                                                pow,
                                                timestep_units='a')

        ######## Run Depletion Calculation
        dep_int.integrate()

        ######## Parse Depletion Calculation Results
        dep_res = openmc.deplete.Results('depletion_results.h5')
        kinf.append(dep_res.get_keff('a')[1].tolist())

        # Calculate and Store Uranium Atom Burnup
        # Accumulate Initial Uranium Atom Concentrations in All Burnable Materials
        U_BOL = 0
        for mat in dep_res.export_to_materials(0):
            if mat.depletable == True:
                U_BOL += dep_res.get_atoms(mat, 'U235', 'atom/cm3', 'a')[1][0]
                U_BOL += dep_res.get_atoms(mat, 'U238', 'atom/cm3', 'a')[1][0]
        
        # Accumulate Final Uranium Atom Concentration in All Burnable Materials
        U_EOL = 0
        for mat in dep_res.export_to_materials(-1):
            if mat.depletable == True:
                U_EOL += dep_res.get_atoms(mat, 'U235', 'atom/cm3', 'a')[1][-1]
                U_EOL += dep_res.get_atoms(mat, 'U238', 'atom/cm3', 'a')[1][-1]
        
        # Calculate Uranium Atom Burnup with Respect to Initial Atom Concentration
        burn.append(((U_BOL - U_EOL) / U_BOL) * 100)

        ######## Delete Depletion Results Container
        del dep_res
    
    # Store Name of Current Moderator
    Mname.append(Modname)

    # Store Moderator to Fuel Density Ratio for Current Moderator
    MF_den_ratio.append(MF_dens)

    # Store Coolant Channel Diameter for Current Moderator
    chnl_dia.append(cc_dia)

    # Store Fuel Element Circumscribing Diameter for Current Moderator
    hex_dia.append(hp_dia)

    # Store Fuel Element Linear Mass Data for Current Moderator
    hex_mass.append(fe_mass)

    # Store Kinf Values for Current Moderator
    K_EVL.append(kinf)
    K_BOL.append([i[0] for i in kinf])
    K_EOL.append([i[-1] for i in kinf])

    # Store Uranium Atom Burnup Values for Current Moderator
    U_burn.append(burn)


############### Save Depletion Calculation Results
with open(Pikname + '.pckl', 'wb') as f:
    pickle.dump(MF_vol, f)
    pickle.dump(dep_sch, f)
    pickle.dump(Mname, f)
    pickle.dump(MF_den_ratio, f)
    pickle.dump(chnl_dia, f)
    pickle.dump(hex_dia, f)
    pickle.dump(hex_mass, f)
    pickle.dump(K_EVL, f)
    pickle.dump(K_BOL, f)
    pickle.dump(K_EOL, f)
    pickle.dump(U_burn, f)


############### Cleanup XML and H5 Files
os.system('rm -rf *.xml')
os.system('rm -rf *.h5')
