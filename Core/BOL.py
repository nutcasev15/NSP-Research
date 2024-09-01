# OpenMC 0.15.0 Python Master Script for Reflector BOL Keff Study


############### Library Imports
import os
import pickle
import numpy as np
import openmc
import openmc.deplete


############## Import Materials Database Symbols and Lookup Functions
from MatDB import *


############### Import Model Builder Function
from BuildModel3D import BuildModel


############### Define Investigation Parameters
# Define Save File Name
Pikname = 'BOL-Ref-0'

# Define Range of Fuel Element Rings in Reactor Core
Rings = np.arange(4, 16 + 1, 1).tolist()

# Define Maximum Fuel Element Linear Power
pow = 275 # W / cm

# Define Radial Reflector Thickness
Ref_thick = 0.1 # cm

# Define Radial Reflector Material ID
Ref_id = 3

# Define List to Store Moderator Name from Model Builder Function
Mname = []

# Define List to Store Reflector Name from Model Builder Function
Rname = []

# Define List to Store Reactor Core Radius from Model Builder Function
core_rad = []

# Define List to Store Reactor Core Height from Model Builder Function
core_hei = []

# Define List to Store Radial Reflector Mass from Model Builder Function
ref_mass = []

# Define List to Store Reactor Core Mass from Model Builder Function
core_mass = []

# Define List to Store BOL Keff Results vs Reflector Thickness
K_BOL = []

# Define List to Store Reactor Core Thermal Power Density
core_tpd = []


############### Import List of Optimum MF Volume Ratios from Kinf BOL Runs
# Tuple Structure: ID, Optimum MF Volume Ratio
ModData = [
    (MATID_BeO,            7.975),
    (MATID_ZrH,            2.092),
    (MATID_YH,             1.546),
    (MATID_MgO_40ZrH,      3.184),
    (MATID_MgO_40YH,       3.008),
    (MATID_BeO_MgO_40ZrH,  3.395),
    (MATID_BeO_MgO_40YH,   3.025),
    (MATID_CaO_CaH,        2.779)
]


############### Run BOL Routine for Each Moderator
for Mod in ModData:
    # Define Temporary Containers for Model Builder and Calculation Outputs
    Modname = []
    Refname = []
    core_r = []
    core_h = []
    ref_m = []
    core_m = []
    k_3D = []
    core_pd = []

    for R in Rings:
        ######## Call Model Builder Function
        Model_Data = BuildModel(Mod[0], Ref_id, Mod[1], R, Ref_thick)

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

        ######## Export and Run Model
        Model_Data[0].export_to_model_xml()
        openmc.run()

        ######## Parse Keff Calculation Results
        k_3D.append(openmc.StatePoint('statepoint.30.h5').keff)

        ######## Calculate and Save Core Power Density for Critical Cores
        core_pd.append((pow * ((3 * R**2) - 3 * R + 1) * Model_Data[4])
                       / (np.pi * Model_Data[3]**2 * Model_Data[4]))

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

    # Store BOL Keff Values for Current Moderator
    K_BOL.append(k_3D)

    # Store Reactor Core Thermal Power Density
    core_tpd.append(core_pd)


############### Save BOL Calculation Results
with open(Pikname + '.pckl', 'wb') as f:
    pickle.dump(Rings, f)
    pickle.dump(pow, f)
    pickle.dump(Ref_thick, f)
    pickle.dump(MF_Mod, f)
    pickle.dump(Mname, f)
    pickle.dump(Rname, f)
    pickle.dump(core_rad, f)
    pickle.dump(core_hei, f)
    pickle.dump(ref_mass, f)
    pickle.dump(core_mass, f)
    pickle.dump(K_BOL, f)
    pickle.dump(core_tpd, f)


############### Cleanup XML and H5 Files
os.system('rm -rf *.xml')
os.system('rm -rf *.h5')
