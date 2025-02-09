# OpenMC 0.15.0 Python Script for BOL MF Ratio Study


############### Library Imports
import os
import pickle
import numpy as np
import openmc


############## Import Materials Database Symbols and Lookup Functions
from MatDB import *


############### Import Model Builder Function
from BuildHexFE2D import BuildModel


############### Define Investigation Parameters
# Define Save File Name
PklName = 'BOL-Ref'

# Define List of 16 MF Volume Ratios between 1 and 10
MF_vol = np.linspace(1, 10, 16).tolist()

# Define List to Store BOL Kinf Results vs MF Volume Ratio
K_BOL = []

# Define List to Store Moderator Name from Model Builder Function
ModName = []

# Define List to Store MF Density Ratio from Model Builder Function
MF_den_ratio = []

# Define List to Store Coolant Channel Diameter from Model Builder Function
channel_dia = []

# Define List to Store Fuel Element Width from Model Builder Function
hex_width = []

# Define List to Store Fuel Element Linear Mass from Model Builder Function
hex_mass = []

# Define List of Moderator Materials to Simulate
ModIDs = [
    MATID_BeO,
    MATID_ZrH,
    MATID_YH,
    MATID_MgO_40ZrH,
    MATID_MgO_40YH,
    MATID_BeO_MgO_40ZrH,
    MATID_BeO_MgO_40YH,
    MATID_CaO_CaH
]

############### Run BOL Routine for Each Moderator
for ID in ModIDs:
    # Define Temporary Containers for Model Builder and Calculation Outputs
    MName = []
    MF_dens = []
    cc_dia = []
    hp_width = []
    fe_mass = []
    kinf = []

    for MF in MF_vol:
        ######## Call Model Builder Function
        Model_Data = BuildModel(ID, MF)

        ######## Save Model Output Data
        MName.append(Model_Data[1])
        MF_dens.append(Model_Data[2])
        cc_dia.append(Model_Data[3])
        hp_width.append(Model_Data[4])
        fe_mass.append(Model_Data[5])

        ######## Clear Run Directory of Previous Results
        os.system('rm -rf *.xml')
        os.system('rm -rf *.h5')

        ######## Run Model and Parse Kinf Results
        kinf.append(openmc.StatePoint(Model_Data[0].run()).keff)

    # Store Name of Current Moderator
    ModName.append(MName)

    # Store Moderator to Fuel Density Ratios for Current Moderator
    MF_den_ratio.append(MF_dens)

    # Store Coolant Channel Diameters for Current Moderator
    channel_dia.append(cc_dia)

    # Store Fuel Element Width for Current Moderator
    hex_width.append(hp_width)

    # Store Fuel Element Linear Mass Data for Current Moderator
    hex_mass.append(fe_mass)

    # Store BOL Kinf Values for Current Moderator
    K_BOL.append(kinf)


############### Save BOL Calculation Results
with open(PklName + '.pkl', 'wb') as f:
    pickle.dump(MF_vol, f)
    pickle.dump(ModName, f)
    pickle.dump(MF_den_ratio, f)
    pickle.dump(channel_dia, f)
    pickle.dump(hex_width, f)
    pickle.dump(hex_mass, f)
    pickle.dump(K_BOL, f)


############### Cleanup XML and H5 Files
os.system('rm -rf *.xml')
os.system('rm -rf *.h5')
