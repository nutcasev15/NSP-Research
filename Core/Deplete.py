# OpenMC 0.13.3 Python Master Script for Reflector Depletion Study


############### Library Imports
import os
import json
import pickle
import numpy as np
import openmc
import openmc.deplete


############### Import Model Builder Function
from OMC_3D import OMC_model


############### Define Investigation Parameters
# Define Save File Name
Pikname = 'Dep-BeO'

# Define Range of Fuel Element Rings in Reactor Core
Rings = np.arange(8, 16 + 1, 1).tolist()

# Define Reactor Core Maximum Thermal Power
pow = 3E6 # W

# Define Radial Reflector Thickness
Ref_thic = 0.1 # cm

# Define Radial Reflector Material ID
Ref_id = 3

# Define Depletion Schedule up to 10 Years
dep_sch = np.logspace(-3, -1, 3) # Setup Schedule to Equilibriate Fission Products
dep_sch = np.append(dep_sch, np.logspace(0, 1, 3)).tolist() # Deplete Until 10 a
dep_steps = np.diff(np.array(dep_sch), prepend=0).tolist() # Get Timesteps for CELI Integrator

# Define Depletion Chain File Name for OpenMC
chain = 'chain_casl_pwr'

# Define List to Store Moderator Name from Model Builder Function
Mname = []

# Define List to Store Reflector Name from Model Builder Function
Rname = []

# Define List to Store Reactor Core Radius from Model Builder Function
core_rad = []

# Define List to Store Reactor Core Height from Model Builder Function
core_hei = []

# Define List to Store Reactor Radial Reflector Mass from Model Builder Function
ref_mass = []

# Define List to Store Reactor Core Mass from Model Builder Function
core_mass = []

# Define Lists to Store Keff Results
K_EVL = []
K_BOL = []
K_EOL = []

# Define List to Store U Atom Burnup
U_burn = []


############### Load Serpent Fission Q Data for Depletion Simulations
with open('serpent_fissq.json', 'r') as f:
    serpent_fission_q = json.load(f)


############### Import List of Optimum MF Volume Ratios from Kinf BOL Runs
MF_Mod = [
    1.000, # Graphite (Dummy Ratio)
    5.016, # ZrH
    3.570, # YH
    1.000, # BeO (Dummy Ratio)
    7.944, # MgO-40YH
    7.845, # BeO.MgO-40YH
    8.676  # CaO-40CaH2
]


############### Run Depletion Routine for Each Moderator
for ID in range(0, 7):
    # Define Temporary Containers for Model Builder and Calculation Outputs
    Modname = []
    Refname = []
    core_r = []
    core_h = []
    ref_m = []
    core_m = []
    k_3D = []
    burn = []

    # Do not Deplete Graphite and BeO Cores
    if ID == 0 or ID == 3:
        continue

    for R in Rings:
        ######## Call Model Builder Function
        Model_Data = OMC_model(ID, Ref_id, MF_Mod[ID], R, Ref_thic)

        ######## Save Model Output Data
        Modname.append(Model_Data[1])
        Refname.append(Model_Data[2])
        core_r.append(Model_Data[3])
        core_h.append(Model_Data[4])
        ref_m.append(Model_Data[5])
        core_m.append(Model_Data[6])

        ######## Clear Run Directory of Previous Results
        os.system('rm -rf *.xml')
        os.system('rm -rf *.h5')

        ######## Override Simulation Settings to Reduce Runtime
        Model_Data[0].settings.batches = 20
        Model_Data[0].settings.inactive = 10
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
        k_3D.append(dep_res.get_keff('a')[1].tolist())

        # Calculate and Store Total Uranium Atom Burnup
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

    # Store Name of Current Reflector
    Rname.append(Refname)

    # Store Reactor Core Radius for Current Moderator
    core_rad.append(core_r)

    # Store Reactor Core Height for Current Moderator
    core_hei.append(core_h)

    # Store Reactor Radial Reflector Mass for Current Moderator
    ref_mass.append(ref_m)

    # Store Reactor Core Mass for Current Moderator
    core_mass.append(core_m)

    # Store Kinf Values for Current Moderator
    K_EVL.append(k_3D)
    K_BOL.append([i[0] for i in k_3D])
    K_EOL.append([i[-1] for i in k_3D])

    # Store Uranium Atom Burnup Values for Current Moderator
    U_burn.append(burn)


############### Save Depletion Calculation Results
with open(Pikname + '.pckl', 'wb') as f:
    pickle.dump(Rings, f)
    pickle.dump(pow, f)
    pickle.dump(Ref_thic, f)
    pickle.dump(dep_sch, f)
    pickle.dump(MF_Mod, f)
    pickle.dump(Mname, f)
    pickle.dump(Rname, f)
    pickle.dump(core_rad, f)
    pickle.dump(core_hei, f)
    pickle.dump(ref_mass, f)
    pickle.dump(core_mass, f)
    pickle.dump(K_EVL, f)
    pickle.dump(K_BOL, f)
    pickle.dump(K_EOL, f)
    pickle.dump(U_burn, f)


############### Cleanup XML and H5 Files
os.system('rm -rf *.xml')
os.system('rm -rf *.h5')
