# OpenMC 0.13.3 Python Master Script for BOL MF Ratio Study


############### Library Imports
import os
import pickle
import numpy as np
import openmc
import openmc.deplete


############### Import Model Builder Function
from OMC_model import OMC_model


############### Define Investigation Parameters
# Define Save File Name
Pikname = 'BOL-Ref'

# Define List of 16 MF Volume Ratios between 1 and 25
MF_vol = np.logspace(0, np.log10(25), 16).tolist()

# Define List to Store BOL Kinf Results vs MF_vol
K_BOL = []

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

############### Run BOL Routine for Each Moderator
for ID in range(0, 7):
    # Define Temporary Containers for Model Builder and Calculation Outputs
    Modname = []
    MF_dens = []
    cc_dia = []
    hp_dia = []
    fe_mass = []
    kinf = []

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

        ######## Export and Run Model
        Model_Data[0].export_to_model_xml()
        openmc.run()

        ######## Parse Kinf Results
        kinf.append(openmc.StatePoint('statepoint.30.h5').keff)

    # Store Name of Current Moderator
    Mname.append(Modname)

    # Store Moderator to Fuel Density Ratios for Current Moderator
    MF_den_ratio.append(MF_dens)

    # Store Coolant Channel Diameters for Current Moderator
    chnl_dia.append(cc_dia)

    # Store Fuel Element Circumscribing Diameters for Current Moderator
    hex_dia.append(hp_dia)

    # Store Fuel Element Linear Mass Data for Current Moderator
    hex_mass.append(fe_mass)

    # Store BOL Kinf Values for Current Moderator
    K_BOL.append(kinf)


############### Save BOL Calculation Results
with open(Pikname + '.pckl', 'wb') as f:
    pickle.dump(MF_vol, f)
    pickle.dump(Mname, f)
    pickle.dump(MF_den_ratio, f)
    pickle.dump(chnl_dia, f)
    pickle.dump(hex_dia, f)
    pickle.dump(hex_mass, f)
    pickle.dump(K_BOL, f)


############### Cleanup XML and H5 Files
os.system('rm -rf *.xml')
os.system('rm -rf *.h5')
