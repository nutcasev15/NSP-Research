# OpenMC 0.15.0 Python Script for MF Ratio Depletion Study


############### Library Imports
import os
import json
import pickle
import numpy as np
import openmc
import openmc.deplete
import openmc.material


############## Import Materials Database Symbols and Lookup Functions
from MatDB import *


############### Import Model Builder Function
from BuildModel2D import BuildModel


############### Import Reactivity Coefficient Calculation Function
from CoefCalc import calculate_reactivity_coef


############### Define Investigation Parameters
# Define Save File Name
PklName = 'Dep-Ref'

# Define Maximum Fuel Element Linear Power
pow = 275 # W / cm


# Define Depletion Schedule upto 45 MWd/kg
deplete_EOL = 45
# Setup Schedule to Equilibrate Fission Products
EQU_sch = np.logspace(-4, -3, 3)
EQU_steps = np.diff(np.array(EQU_sch), prepend=0).tolist()

# Setup Schedule to Deplete until MOL
MOL_sch = np.logspace(-2, np.log10(deplete_EOL / 2), 3)
MOL_steps = np.diff(np.array(MOL_sch), prepend=EQU_sch[-1]).tolist()

# Setup Schedule to Deplete until EOL
EOL_sch = np.logspace(np.log10(deplete_EOL * 0.667), np.log10(deplete_EOL), 3)
EOL_steps = np.diff(np.array(EOL_sch), prepend=MOL_sch[-1]).tolist()

# Store Complete Depletion Schedule for Post Processing
dep_sch = EQU_sch.tolist() + MOL_sch.tolist() + EOL_sch.tolist()

# Define Depletion Chain File Name for OpenMC
chain = 'chain_casl_pwr'

# Define List to Store Moderator Name from Model Builder Function
Mname = []

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

# Define Lists to Store Reactivity Calculation Results
Coef_BOL = []
Coef_XST = []
Coef_MOL = []
Coef_EOL = []


############### Define List of Moderator Materials to Simulate
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

############### Import Serpent Fission Q Data for Depletion Simulations
with open('serpent_fission_q.json', 'r') as f:
    serpent_fission_q = json.load(f)


############### Define Suite of Reactivity Coefficient Calculations
def find_model_coef(model : openmc.Model,
                    mod : openmc.Material, fuel : openmc.Material, cool : openmc.Material,
                    percent_change = 10.0):
    # Define Material Density Updater Function
    def model_density_updater(model : openmc.Model, deviation, mat : openmc.Material):
        dens = (1.0 + deviation) * get_materials_by_name(mat.name)[0].get_mass_density()
        model.update_densities([mat.id], dens, density_units='g/cm3')
        return dens

    # Define Cell Temperature Updater Function
    def model_temperature_updater(model : openmc.Model, deviation, mat : openmc.Material):
        temp = (1.0 + deviation) * mat.temperature
        model.update_cell_temperatures([mat.name], temp)

        if mat.name is get_material_by_id(MATID_UN).name:
            model.update_densities([mat.id], get_UN_Ideal_Dens(temp), density_units='g/cm3')
        elif mat.name is get_material_by_id(MATID_HeXe).name:
            model.update_densities([mat.id], get_HeXe_40MM_2MPa_Dens(temp), density_units='g/cm3')

        return temp

    # Define Percent Deviation Input Values
    input_dev = np.arange(start=-percent_change, stop=(percent_change + 1.0), step=percent_change)
    input_dev /= 100
    input_dev = input_dev.tolist()

    # Calculate Moderator Reactivity Coefficients
    original_mod_dens = mod.get_mass_density()
    original_mod_temp = mod.temperature
    mod_temp_values, mod_temp_coef = calculate_reactivity_coef(model=model,
                                                               param_list=input_dev,
                                                               model_updater=model_temperature_updater,
                                                               updater_static_args={'mat' : mod},
                                                               print_output=False)
    model.update_cell_temperatures([mod.name], original_mod_temp)
    mod_dens_values, mod_dens_coef = calculate_reactivity_coef(model=model,
                                                               param_list=input_dev,
                                                               model_updater=model_density_updater,
                                                               updater_static_args={'mat' : mod},
                                                               print_output=False)
    model.update_densities([mod.id], original_mod_dens, density_units='g/cm3')


    # Calculate Fuel Temperature Reactivity Coefficients
    original_fuel_temp = fuel.temperature
    fuel_temp_values, fuel_temp_coef = calculate_reactivity_coef(model=model,
                                                               param_list=input_dev,
                                                               model_updater=model_temperature_updater,
                                                               updater_static_args={'mat' : fuel},
                                                               print_output=False)
    model.update_cell_temperatures([fuel.name], original_fuel_temp)
    model.update_densities([fuel.id], get_UN_Ideal_Dens(original_fuel_temp), density_units='g/cm3')

    # Calculate HeXe Coolant Temperature Reactivity Coefficients
    original_cool_temp = cool.temperature
    cool_temp_values, cool_temp_coef = calculate_reactivity_coef(model=model,
                                                               param_list=input_dev,
                                                               model_updater=model_temperature_updater,
                                                               updater_static_args={'mat' : cool},
                                                               print_output=False)
    model.update_cell_temperatures([cool.name], original_cool_temp)
    model.update_densities([cool.id], get_HeXe_40MM_2MPa_Dens(original_cool_temp), density_units='g/cm3')

    return [(mod_dens_values, mod_dens_coef), (mod_temp_values, mod_temp_coef), \
        (fuel_temp_values, fuel_temp_coef), (cool_temp_values, cool_temp_coef)]


############### Run Depletion Routine for Each Moderator
for Mod in ModData:
    ######## Call Model Builder Function
    Model_Data = BuildModel(Mod[0], Mod[1])

    ######## Clear Run Directory of Previous Results
    os.system('rm -rf *.xml')
    os.system('rm -rf *.h5')

    ####### Override Simulation Settings to Reduce Runtime
    Model_Data[0].settings.batches = 15
    Model_Data[0].settings.inactive = 5
    Model_Data[0].settings.particles = 100

    Model_Data[0].init_lib()

    BOL_Coef = find_model_coef(Model_Data[0],
                               Model_Data[6],
                               Model_Data[7],
                               Model_Data[8],
                               30.0)

    Model_Data[0].deplete(EQU_steps,
                          method='celi',
                          operator_kwargs={'chain_file' : chain,
                                           'fission_yield_mode' : 'average',
                                           'fission_q' : serpent_fission_q,
                                           'diff_burnable_mats': True},
                          timestep_units='MWd/kg',
                          power=pow)

    # Accumulate Initial Uranium Atom Concentrations in All Burnable Materials
    U_BOL = 0
    dep_res = openmc.deplete.Results('depletion_results.h5')
    for mat in dep_res.export_to_materials(burnup_index=0):
        if mat.depletable == True:
            U_BOL += dep_res.get_atoms(mat, 'U235', 'atom/cm3', 'a')[1][0]
            U_BOL += dep_res.get_atoms(mat, 'U238', 'atom/cm3', 'a')[1][0]

    XST_Coef = find_model_coef(Model_Data[0],
                               Model_Data[6],
                               Model_Data[7],
                               Model_Data[8],
                               30.0)

    Model_Data[0].deplete(MOL_steps,
                          method='celi',
                          operator_kwargs={'chain_file' : chain,
                                           'fission_yield_mode' : 'average',
                                           'fission_q' : serpent_fission_q,
                                           'diff_burnable_mats': True},
                          timestep_units='MWd/kg',
                          power=pow)

    MOL_Coef = find_model_coef(Model_Data[0],
                            Model_Data[6],
                            Model_Data[7],
                            Model_Data[8],
                            30.0)

    Model_Data[0].deplete(EOL_steps,
                          method='celi',
                          operator_kwargs={'chain_file' : chain,
                                           'fission_yield_mode' : 'average',
                                           'fission_q' : serpent_fission_q,
                                           'diff_burnable_mats': True},
                          timestep_units='MWd/kg',
                          power=pow)

    EOL_Coef = find_model_coef(Model_Data[0],
                               Model_Data[6],
                               Model_Data[7],
                               Model_Data[8],
                               30.0)

    Model_Data[0].finalize_lib()

    ######## Parse Depletion Calculation Results
    dep_res = openmc.deplete.Results('depletion_results.h5')
    kinf = dep_res.get_keff('a')[1].tolist()

    # Calculate and Store Uranium Atom Burnup
    # Accumulate Final Uranium Atom Concentration in All Burnable Materials
    U_EOL = 0
    for mat in dep_res.export_to_materials(burnup_index=-1):
        if mat.depletable == True:
            U_EOL += dep_res.get_atoms(mat, 'U235', 'atom/cm3', 'a')[1][-1]
            U_EOL += dep_res.get_atoms(mat, 'U238', 'atom/cm3', 'a')[1][-1]

    # Calculate Uranium Atom Burnup with Respect to Initial Atom Concentration
    burn = ((U_BOL - U_EOL) / U_BOL) * 100

    ######## Delete Depletion Results Container
    del dep_res

    # Store Name of Current Moderator
    Mname.append(Model_Data[1])

    # Store Moderator to Fuel Density Ratio for Current Moderator
    MF_den_ratio.append(Model_Data[2])

    # Store Coolant Channel Diameter for Current Moderator
    channel_dia.append(Model_Data[3])

    # Store Fuel Element Circumscribing Diameter for Current Moderator
    hex_dia.append(Model_Data[4])

    # Store Fuel Element Linear Mass Data for Current Moderator
    hex_mass.append(Model_Data[5])

    # Store Kinf Values for Current Moderator
    K_EVL.append(kinf)
    K_BOL.append([i[0] for i in kinf])
    K_EOL.append([i[-1] for i in kinf])

    # Store Uranium Atom Burnup Values for Current Moderator
    U_burn.append(burn)

    Coef_BOL.append(BOL_Coef)
    Coef_XST.append(XST_Coef)
    Coef_MOL.append(MOL_Coef)
    Coef_EOL.append(EOL_Coef)

    del Model_Data


############### Save Depletion Calculation Results
with open(PklName + '.pkl', 'wb') as f:
    pickle.dump(ModData, f)
    pickle.dump(dep_sch, f)
    pickle.dump(Mname, f)
    pickle.dump(MF_den_ratio, f)
    pickle.dump(channel_dia, f)
    pickle.dump(hex_dia, f)
    pickle.dump(hex_mass, f)
    pickle.dump(K_EVL, f)
    pickle.dump(K_BOL, f)
    pickle.dump(K_EOL, f)
    pickle.dump(U_burn, f)
    pickle.dump(Coef_BOL, f)
    pickle.dump(Coef_XST, f)
    pickle.dump(Coef_MOL, f)
    pickle.dump(Coef_EOL, f)


############### Cleanup XML and H5 Files
os.system('rm -rf *.xml')
os.system('rm -rf *.h5')
