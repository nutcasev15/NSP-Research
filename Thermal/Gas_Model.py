# Bounding Fuel Element Linear Power Calculator for Gas Cooled NSP Systems
# Powered by NumPy 1.26.4 and SciPy 1.14.0
#
# Refer ntrs.nasa.gov/citations/20060056311 for HeXe Coolant Data
# Refer doi.org/10.1016/0009-2614(75)80286-7 for He Coolant Data
# Refer iaea.org/publications/7965/thermophysical-properties-of\
#       -materials-for-nuclear-engineering-a-tutorial-and-collection-of-data
#       for UN Integrated Thermal Conductivity Data
# Refer osti.gov/biblio/532973 for MA956 ODS Steel Thermal Conductivity Data


############### Library Imports
from os import system, name
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import root_scalar

import itertools
import matplotlib.pyplot as plt
import matplotlib.markers as mrk
import matplotlib.ticker as tck


############### Set Script Inputs
# UN Fuel Pellet Radius
Rf = 7.9 # mm

# MA 956 ODS Steel Radius
Rc = 8.4 # mm

# Maximum Fuel Centreline Temperature
Tm = 1600 # K

# Bulk Coolant Temperature Range
Tb = np.arange(900, 1500, 100) # K

# Coolant Mass Flux Range
G = np.arange(100, 1100, 25) # kg / (m^2 * s)

# Initialise Fuel Linear Power Rating
Q = np.NaN * np.ones((Tb.shape[0], G.shape[0])) # W / m


############### Import and Define Material Thermal Data
# Coolant HTC Values from Johnson-2006 and Taylor-1988
HTCData = np.genfromtxt('HeXe_HTC.csv', delimiter=',')

# Coolant HTC Values from Jain-1975 and Taylor-1988
# HTCData = np.genfromtxt('He_HTC.csv', delimiter=',')

# Replace NaN Values with a Placeholder for Interpolation
HTCData[np.isnan(HTCData)] = 1

# Setup Interpolator for HTC Values
get_HTC = RegularGridInterpolator((HTCData[1:, 0],  # T in K
                                   HTCData[0, 1:]), # G in kg / (m^2 * s)
                                   HTCData[1:, 1:], # HTC in W / (m^2 * K)
                                   'cubic')

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
def get_TH_Ts(x, tb, g):
    return (tb + x * ((np.log(Rc / Rf) / (2 * np.pi * get_k_MA956(tb))) \
            + (1 / (2 * np.pi * Rc * 1E-3 * get_HTC((tb, g))))))


############### Calculate Linear Power Values
Tol = 1E-12
Q_Low = 25E2 # W / m
T_It = np.nditer(Tb, flags=['f_index'])
G_It = np.nditer(G, flags=['f_index'])

for tb in T_It:
    for g in G_It:
        # Find Fuel Pin Linear Power Rating in W / m
        # Raise Function's Output to 1.39 to Avoid Imaginary Values
        def opt_func(x):
            return (np.real(get_UN_Ts(x)**1.39)
                    - (get_TH_Ts(x, tb, g)**1.39))

        try:
            res = root_scalar(opt_func, bracket=[100, 1E6], xtol=Tol,
                            method='toms748')
        except ValueError:
            continue

        # Handle Case for No Convergence by Rejecting the Result
        # Also Reject Solutions with Very Low Linear Powers
        if res.converged is not True or res.root < Q_Low:
            continue

        # Save the Root Finder Results as the Calculated Linear Power
        Q[T_It.index, G_It.index] = res.root

    # Reset Mass Flux Values Iterator for Next Temperature Value
    G_It.reset()

    # Clear Console Output and Show Solution Progress
    if name == 'nt':
        system('cls')
    else:
        system('clear')

    print('Linear Power Calculation Progress >> {:2.1f} %'.format(
            ((100 * T_It.index) / Tb.shape[0])))


############### Plot Solutions
####### Define Set of Markers for Plotting
marks = itertools.cycle(mrk.MarkerStyle.filled_markers)

###### Define Set of Line Styles for Plotting
lines = itertools.cycle(['-', '--', ':', '-.'])

####### Reactor Linear Power with Fixed Range for Comparison
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title('Gas-Cooled Reactor Maximum UN Pellet Linear Power')
ax.set_xlabel('Core Mass Flux (kg/m$^2$s)')
ax.grid()

for i in range(0, Tb.shape[0]):
    # Plot Calculated Linear Power Data for Coolant Temperature
    ax.plot(G, Q[i, :], label='HeXe at {:.0f} K'.format(Tb[i]),
            marker=next(marks), linestyle=next(lines))

# Update Y Axis Values to W / cm
ax.yaxis.set_major_locator(tck.MultipleLocator(50E2))
ax.yaxis.set_minor_locator(tck.MultipleLocator(25E2))
ax.yaxis.set_major_formatter(lambda x, _: str(int(x * 1E-2)))
ax.set_ylim(bottom=Q_Low)
ax.set_ylabel('Reactor Maximum Linear Power (W/cm)')

# Finalise and Save Reactor Linear Power Plot
ax.legend(loc='upper left', prop={'size' : 8})
fig.tight_layout()
fig.savefig(f'{Rf:2.0f}' + '_mm_Gas_Linear_Power.pdf', format='pdf')

# Clear Figure and Axes
fig.clear()
ax.clear()
