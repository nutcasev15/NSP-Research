# OpenMC 0.15.2 Python Script for BOL MF Ratio Study


############### Library Imports
import os
from sys import exit
import pickle
import numpy as np
import openmc


############## Import Materials Database Symbols and Lookup Functions
from MatDB import *


############### Import Model Builder Function
from BuildHexFE2D import BuildModel


############### Define Investigation Parameters
# Define Save File Name
pkl_name = 'BOL-Ref'

# Define List of 16 MF Volume Ratios between 1 and 12
MF_vol = np.logspace(np.log10(1), np.log10(12), 16).tolist()

# Define List to Store BOL Kinf Results vs MF Volume Ratio
K_BOL = []

# Define List to Store Moderator Name from Model Builder Function
mod_name = []

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
  MATID_Graphite,
  MATID_ZrH,
  MATID_YH,
  MATID_BeO,
  MATID_MgO_40YH,
  MATID_BeO_MgO_40YH,
  MATID_CaO_CaH
]

############### Run BOL Routine for Each Moderator
for ID in ModIDs:
  # Define Temporary Containers for Model Builder and Calculation Outputs
  name = []
  MF_dens = []
  cc_dia = []
  hp_width = []
  fe_mass = []
  kinf = []

  for MF in MF_vol:
    try:
      ######## Call Model Builder Function
      print('>> Building Model <<')
      Model_Data = BuildModel(MID=ID, MFR=MF)

      ######## Save Model Output Data
      name.append(Model_Data['Moderator Material'].name)
      MF_dens.append(Model_Data['Density Ratio'])
      cc_dia.append(Model_Data['Channel Diameter'])
      hp_width.append(Model_Data['FE Width'])
      fe_mass.append(Model_Data['FE Mass'])

      ######## Clear Run Directory of Previous Results
      print('>> Clearing Previous Result Files <<')
      os.system('rm -rf *.xml')
      os.system('rm -rf *.h5')

      ######## Run Model and Parse Kinf Results
      print('>> Moderator: {:s} with Volume Ratio: {:.3f} <<'.format(
        name[-1],
        MF
      ))
      result = openmc.StatePoint(Model_Data['Model'].run(output=False))
      kinf.append(
        result.keff
      )

      ######## Output Quick Diagnostics
      print('>> Multiplication Factor: {:1.6f} +/- {:1.6} <<'.format(
        result.keff.n,
        result.keff.s
      ))
      print('>> Total Simulation Time: {:3.1f} s <<'.format(
        result.runtime['total']
      ))
      print('\n')
    except KeyboardInterrupt:
      ######## Abort Simulation without Saving Results
      exit()
    except:
      ######## Stop Simulation on Error Preserving Existing Results
      break

  # Store Name of Current Moderator
  mod_name.append(name)

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
with open(pkl_name + '.pkl', 'wb') as f:
    pickle.dump(MF_vol, f)
    pickle.dump(mod_name, f)
    pickle.dump(MF_den_ratio, f)
    pickle.dump(channel_dia, f)
    pickle.dump(hex_width, f)
    pickle.dump(hex_mass, f)
    pickle.dump(K_BOL, f)


############### Cleanup XML and H5 Files
os.system('rm -rf *.xml')
os.system('rm -rf *.h5')
