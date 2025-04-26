# OpenMC 0.15.2 Python Script for MF Ratio Depletion Study


############### Library Imports
import os
from sys import exit
import json
import pickle
import numpy as np
import openmc
from openmc.deplete import Results


############## Import Materials Database Symbols and Lookup Functions
from MatDB import *


############### Import Model Builder Function
from BuildHexFE2D import BuildModel


############### Import Reactivity Coefficient Calculation Function
from CoefCalc import calculate_reactivity_coef


############### Define Investigation Parameters
# Define Save File Name
pkl_name = 'Dep-Ref'

# Define Maximum Fuel Element Linear Power
pow = 250 # W / cm

# Define Temperature Deviation in Percent for Reactivity Calculation
dev = 15.0

# Store Cumulative Depletion Schedule in Days for Post Processing
dep_sch = [
  0.3,         # 8 Hours
  0.6,
  1.0,         # Xenon Stabilisation Time
  3.0,
  10.0,
  30.0,        # 1 Month
  90.0,
  180.0,
  365.0,       # 1 Year
  5.0 * 365.0, # Middle of Life
  10.0 * 365.0 # EOL at 10 Years
]

# Convert Depletion Schedule to Timesteps for Solver
dep_steps = np.diff(dep_sch).tolist()

# Define Depletion Chain File Name for OpenMC
chain = 'chain_casl_pwr'

# Define List to Store Moderator Name from Model Builder Function
mod_name = []

# Define List to Store MF Density Ratio from Model Builder Function
MF_den_ratio = []

# Define List to Store Coolant Channel Diameter from Model Builder Function
channel_dia = []

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

# # Define Lists to Store Reactivity Calculation Results
coef_BOL = []
coef_EOL = []


############### Define List of Moderator Materials to Simulate
# Dictionary Structure: ID, Optimum MF Volume Ratio
mod_data = [
    {'ID': MATID_Graphite,       'MFR' : 12.00},
    {'ID': MATID_ZrH,            'MFR' : 1.775},
    {'ID': MATID_YH,             'MFR' : 1.247},
    {'ID': MATID_BeO,            'MFR' : 12.00},
    {'ID': MATID_MgO_40YH,       'MFR' : 2.603},
    {'ID': MATID_BeO_MgO_40YH,   'MFR' : 2.708},
    {'ID': MATID_CaO_CaH,        'MFR' : 2.321}
]

############### Import Serpent Fission Q Data for Depletion Simulations
with open('serpent_fission_q.json', 'r') as f:
    serpent_fission_q = json.load(f)


############### Define Suite of Reactivity Coefficient Calculations
def find_model_coef(builder_output : dict, percent_change : float = 10.0):
  # Define Temperature Updater Function
  def model_updater(model : openmc.Model, deviation,
                    MID : int, MFR : float, TIT : float):
    # Modify Temperatures and Rebuild Model
    outlet_temp = (1.0 + deviation) * TIT
    updated_model = BuildModel(MID=MID, MFR=MFR, TIT=outlet_temp)

    # Update Cell Temperatures
    # Modify Only Cells With Calculated Cell Temperatures
    for cell in updated_model['Model']._cells_by_id.values():
      cell : openmc.Cell
      if cell.name != '' and isinstance(cell.temperature, float):
        model.update_cell_temperatures([cell.name], cell.temperature)

    # Update Material Nuclide Densities
    # Nuclide Ratios and Composition are Not Modified
    for mat in updated_model['Model'].materials:
      mat : openmc.Material
      if mat.name != '':
        model.update_densities([mat.name], mat.density, 'g/cm3')

    return outlet_temp

  # Define Percent Deviation Input Values
  input_dev = np.arange(
    start=-percent_change,
    stop=(percent_change + 1.0),
    step=percent_change
  )
  input_dev /= 100
  input_dev = input_dev.tolist()

  # Run Reactivity Coefficient Calculation in Memory
  model_coef = calculate_reactivity_coef(
    model=builder_output['Model'],
    param_list=input_dev,
    model_updater=model_updater,
    updater_static_args={
      'MID' : builder_output['Moderator ID'],
      'MFR' : builder_output['Volume Ratio'],
      'TIT' : builder_output['Outlet Temperature']
    },
    started_initialized=True,
    print_output=False
  )

  # Reset Temperatures and Densities in Memory
  # All Other Model Data Left Unchanged from Input State
  model_updater(builder_output['Model'],
    deviation=0.0,
    MID=builder_output['Moderator ID'],
    MFR=builder_output['Volume Ratio'],
    TIT=builder_output['Outlet Temperature']
  )

  return model_coef


############### Run Depletion Routine for Each Moderator
for mod in mod_data:
  try:
    ######## Call Model Builder Function
    print('>> Building Model <<')
    model_data = BuildModel(MID=mod['ID'], MFR=mod['MFR'])
    model : openmc.Model = model_data['Model']
    print('>> Moderator: {:s} with Volume Ratio: {:.3f} <<'.format(
      model_data['Moderator Material'].name,
      mod['MFR']
    ))

    ######## Clear Run Directory of Previous Results
    print('>> Clearing Previous Result Files <<')
    os.system('rm -rf *.xml')
    os.system('rm -rf *.h5')

    ######## Initialise Model in Memory
    print('>> Setting Up Simulation in Memory <<')
    model.init_lib(output=False)

    ######## Calculate Temperature Reactivity Coefficients at BOL
    print('>> Calculating Reactivity Coefficients at BOL <<')
    BOL_coef = find_model_coef(model_data, dev)

    ######## Deplete Model According to Schedule in Memory
    print('>> Depleting Model in Memory <<')
    model.deplete(
      dep_steps,
      method='celi',
      operator_kwargs={
        'chain_file' : chain,
        'fission_yield_mode' : 'average',
        'fission_q' : serpent_fission_q,
        'diff_burnable_mats': True
      },
      timestep_units='d',
      power=pow,
      output=False
    )

    ######## Calculate Temperature Reactivity Coefficients after Depletion
    print('>> Calculating Reactivity Coefficients at EOL <<')
    EOL_coef = find_model_coef(model_data, dev)

    ######## Stop Simulation and Deallocate Memory
    print('>> Finalising Simulation and Deallocating Memory <<')
    print('>> Moderator: {:s} with Volume Ratio: {:.3f} <<'.format(
      model_data['Moderator Material'].name,
      mod['MFR']
    ))
    model.finalize_lib()
    openmc.reset_auto_ids()

    ######## Parse Depletion Calculation Results
    dep_res = Results('depletion_results.h5')
    kinf = dep_res.get_keff('d')[1].tolist()

    # Calculate and Store Uranium Atom Burnup
    # Accumulate Final Uranium Atom Concentration in All Burnable Materials
    U_BOL = 0
    U_EOL = 0
    for mat in dep_res.export_to_materials(-1):
      if mat.depletable == True:
        U_BOL += dep_res.get_atoms(mat, 'U235', 'atom/cm3', 'd')[1][0]
        U_BOL += dep_res.get_atoms(mat, 'U238', 'atom/cm3', 'd')[1][0]
        U_EOL += dep_res.get_atoms(mat, 'U235', 'atom/cm3', 'd')[1][-1]
        U_EOL += dep_res.get_atoms(mat, 'U238', 'atom/cm3', 'd')[1][-1]

    # Calculate Uranium Atom Burnup with Respect to Initial Atom Concentration
    burn = ((U_BOL - U_EOL) / U_BOL) * 100
    print('>> Uranium Burnup: {:1.3f} % <<'.format(burn))
    print('\n')

    # Delete Depletion Results Container
    del dep_res
  except KeyboardInterrupt:
    ######## Abort Simulation without Saving Results
    exit()
  except:
    ######## Stop Simulation on Error Preserving Existing Results
    break

  # Store Name of Current Moderator
  mod_name.append(model_data['Moderator Material'].name)

  # Store Moderator to Fuel Density Ratio for Current Moderator
  MF_den_ratio.append(model_data['Density Ratio'])

  # Store Coolant Channel Diameter for Current Moderator
  channel_dia.append(model_data['Channel Diameter'])

  # Store Fuel Element Circumscribing Diameter for Current Moderator
  hex_dia.append(model_data['FE Width'])

  # Store Fuel Element Linear Mass Data for Current Moderator
  hex_mass.append(model_data['FE Mass'])

  # Store Kinf Values for Current Moderator
  K_EVL.append(kinf)
  K_BOL.append(kinf[0])
  K_EOL.append(kinf[-1])

  # Store Uranium Atom Burnup Values for Current Moderator
  U_burn.append(burn)

  # Store Temperature Reactivity Coefficients for Current Moderator
  coef_BOL.append(BOL_coef)
  coef_EOL.append(EOL_coef)

  # Delete Python Instance of Simulation Data
  del model_data


############### Save Depletion Calculation Results
with open(pkl_name + '.pkl', 'wb') as f:
    pickle.dump(mod_data, f)
    pickle.dump(dep_sch, f)
    pickle.dump(mod_name, f)
    pickle.dump(MF_den_ratio, f)
    pickle.dump(channel_dia, f)
    pickle.dump(hex_dia, f)
    pickle.dump(hex_mass, f)
    pickle.dump(K_EVL, f)
    pickle.dump(K_BOL, f)
    pickle.dump(K_EOL, f)
    pickle.dump(U_burn, f)
    pickle.dump(coef_BOL, f)
    pickle.dump(coef_EOL, f)


############### Cleanup XML and H5 Files
os.system('rm -rf *.xml')
os.system('rm -rf *.h5')
