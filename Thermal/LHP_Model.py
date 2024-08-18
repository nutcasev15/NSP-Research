# Bounding Fuel Element Linear Power Calculator for Sodium Cooled NSP Systems
# Powered by NumPy 1.26.4 and SciPy 1.14.0
#
# Refer iaea.org/publications/7965/thermophysical-properties-of\
#       -materials-for-nuclear-engineering-a-tutorial-and-collection-of-data
#       for UN Integrated Thermal Conductivity Data
# Refer osti.gov/biblio/532973 for MA956 ODS Steel Thermal Conductivity Data


############### Library Imports
from os import system, name
import numpy as np
from scipy.optimize import root_scalar

import matplotlib.pyplot as plt
import matplotlib.ticker as tck


############### Set Script Inputs
# UN Fuel Pellet Radius
Rf = 7.9 # mm

# MA 956 ODS Steel Radius
Rc = 8.4 # mm

# Maximum Fuel Centreline Temperature
Tm = 1600 # K

# Sodium Saturated Vapour Temperature Range
Tb = np.arange(800, 1600, 25) # K

# Initialise Fuel Linear Power Rating
Q = np.NaN * np.ones(Tb.shape[0]) # W / m


############### Import and Define Material Thermal Data
# Define Function for MA956 ODS Steel Thermal Conductivity
def get_k_MA956(T):
    return ((1 / 60) * (T - 273) + 10)


############### Define Optimisation Functions
# Define Helper Function to Solve for Fuel Pellet Surface Temperature
# Derived by Integrating the Equation for UN Thermal Conductivity
def get_UN_Ts(x):
    return ((Tm**1.39 - (1.39 / 1.41) * (x / (4 * np.pi)))**(1 / 1.39))

# Define Helper Function to Solve for Fuel Pellet Surface Temperature
# Derived from Radial Heat Conduction and Convection Equations
def get_TH_Ts(x, tb):
    return (tb + x * ((np.log(Rc / Rf) / (2 * np.pi * get_k_MA956(tb)))))


############### Calculate Linear Power Values
Tol = 1E-12
T_It = np.nditer(Tb, flags=['f_index'])

for tb in T_It:
    # Find Fuel Pin Linear Power Rating in W / m
    # Raise Function's Output to 1.39 to Avoid Imaginary Values
    def opt_func(x):
        return (np.real(get_UN_Ts(x)**1.39) - (get_TH_Ts(x, tb)**1.39))

    try:
        res = root_scalar(opt_func, bracket=[100, 1E6], xtol=Tol,
                        method='toms748')
    except ValueError:
        continue

    # Handle Case for No Convergence by Rejecting the Result
    if res.converged is not True:
        continue

    # Save the Root Finder Results as the Calculated Linear Power
    Q[T_It.index] = res.root

    # Clear Console Output and Show Solution Progress
    if name == 'nt':
        system('cls')
    else:
        system('clear')

    print('Linear Power Calculation Progress >> {:2.1f} %'.format(
            ((100 * T_It.index) / Tb.shape[0])))


############### Plot Solutions
####### Reactor Linear Power with Fixed Range for Comparison
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title('Sodium-Cooled Reactor Maximum UN Pellet Linear Power')
ax.set_xlabel('Sodium Vapour Outlet Temperature (K)')
ax.grid()

# Plot Calculated Linear Power Data for Coolant Temperature
ax.plot(Tb, Q)

# Update Y Axis Values to W / cm
ax.yaxis.set_major_locator(tck.MultipleLocator(200E2))
ax.yaxis.set_minor_locator(tck.MultipleLocator(100E2))
ax.yaxis.set_major_formatter(lambda x, _: str(int(x * 1E-2)))
ax.set_ylim(bottom=50E2)
ax.set_ylabel('Reactor Maximum Linear Power (W/cm)')

# Finalise and Save Reactor Linear Power Plot
fig.tight_layout()
fig.savefig(f'{Rf:2.0f}' + '_mm_LHP_Linear_Power.pdf', format='pdf')

# Clear Figure and Axes
fig.clear()
ax.clear()
