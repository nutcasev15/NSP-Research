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
from numpy.polynomial.polynomial import Polynomial
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
Rc = 8.9 # mm

# Target Core Coolant Mass Flux from G_lim.py
G_tar = 600 # kg / (m^2 * s)

# Target Core Linear Power from G_lim.py
Q_tar = 250E2 # W / m

# Maximum Fuel Centreline Temperature
Tm = 1600 # K

# Bulk Coolant Temperature Range
Tb = np.arange(800, 1500, 100) # K

# Coolant Mass Flux Range
G = np.arange(300, 1100, 25) # kg / (m^2 * s)

# Initialise Fuel Linear Power Rating
Q = np.NaN * np.ones((Tb.shape[0], G.shape[0])) # W / m

# Initialise Expected Temperatures from G_lim.py Inputs
# Clad Surface
CladTemp = np.NaN * np.ones((Tb.shape[0], )) # K
# Fuel Pellet Surface
FuelTemp = np.NaN * np.ones((Tb.shape[0], )) # K
# Fuel Pellet Centreline
CentreTemp = np.NaN * np.ones((Tb.shape[0], )) # K


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
# Define Helper Function to Solve for Fuel Pellet Centreline Temperature
# Derived by Integrating the Equation for UN Thermal Conductivity
def get_UN_Tm(x, ts):
    return ((ts**1.39 + (1.39 / 1.41) * (x / (4 * np.pi)))**(1 / 1.39))

# Define Helper Function to Solve for Fuel Pellet Surface Temperature
# Derived by Integrating the Equation for UN Thermal Conductivity
def get_UN_Ts(x):
    return ((Tm**1.39 - (1.39 / 1.41) * (x / (4 * np.pi)))**(1 / 1.39))

# Define Helper Function to Solve for Fuel Pellet Surface Temperature
# Derived from Radial Heat Conduction and Convection Equations
def get_TH_Ts(x, tb, g):
    return (tb + x * ((np.log(Rc / Rf) / (2 * np.pi * get_k_MA956(tb))) \
            + (1 / (2 * np.pi * Rc * 1E-3 * get_HTC((tb, g))))))

# Define Helper Function to Solve for Clad Surface Temperature
# Derived from Radial Heat Convection Equations
def get_TH_Tc(x, tb, g):
    return (tb + x * (1 / (2 * np.pi * Rc * 1E-3 * get_HTC((tb, g)))))


############### Calculate Linear Power Values
Tol = 1E-12
Q_Low = 25E2 # W / m
T_It = np.nditer(Tb, flags=['f_index'])
G_It = np.nditer(G, flags=['f_index'])

for tb in T_It:
    # Calculate Temperatures at Expected Mass Flux and Linear Power
    CladTemp[T_It.index] = get_TH_Tc(Q_tar, tb, G_tar)
    FuelTemp[T_It.index] = get_TH_Ts(Q_tar, tb, G_tar)
    CentreTemp[T_It.index] = get_UN_Tm(Q_tar, FuelTemp[T_It.index])

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


############### Output Fuel Element Temperatures as Curve Fits
# Curve Fit to 4th Order Polynomial as a Function of Coolant Temperature
CentreTempCoef = Polynomial.fit(Tb, CentreTemp, 4).convert().coef
FuelTempCoef = Polynomial.fit(Tb, FuelTemp, 4).convert().coef
CladTempCoef = Polynomial.fit(Tb, CladTemp, 4).convert().coef

# Output Equations and Validity Criteria for Curve Fits
print('\n###############')
print('From G_lim.py')
print('{:10} {:10.3f} kg/m^2*s'.format('Target G:', G_tar))
print('{:10} {:10.3f} W/m'.format('Target q\':', Q_tar))
print('\nFuel Element Fitted Temperatures vs Bulk Coolant Temperature')
Equation = '{:2} = {:1.4E} * Tb**4 ' + \
    '+ {:1.4E} * Tb**3 ' + \
    '+ {:1.4E} * Tb**2 ' + \
    '+ {:1.4E} * Tb ' + \
    '+ {:1.4E}'
print(Equation.format(
    'Tm',
    CentreTempCoef[4],
    CentreTempCoef[3],
    CentreTempCoef[2],
    CentreTempCoef[1],
    CentreTempCoef[0]
))
print(Equation.format(
    'Ts',
    FuelTempCoef[4],
    FuelTempCoef[3],
    FuelTempCoef[2],
    FuelTempCoef[1],
    FuelTempCoef[0]
))
print(Equation.format(
    'Tc',
    CladTempCoef[4],
    CladTempCoef[3],
    CladTempCoef[2],
    CladTempCoef[1],
    CladTempCoef[0]
))


############### Plot Solution
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
print('\n###############')
print('Reactor Linear Power Plot Saved as PDF')

# Clear Figure and Axes
fig.clear()
ax.clear()
