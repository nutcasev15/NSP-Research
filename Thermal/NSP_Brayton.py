# Bounding Thermodynamic Parameter Calculator for Simple Gas-Cooled NSP Systems
# Powered by NumPy 1.26.4 and SciPy 1.14.0
#
# Refer NSP_Brayton.pdf for Equations
# Refer ec.europa.eu/research/participants/documents/downloadPublic?\
#       documentIds=080166e5c802fe77&appId=PPGMS for Radiator Parameters
# Refer doi.org/10.2514/1.27664 for He-Xe Coolant Data


############### Library Imports
from os import system, name
import numpy as np
from scipy.optimize import root_scalar

import matplotlib.pyplot as plt


############### Define Natural Constants
# Import Physical Constants in SI Units
R = 8.314462618
s_b = 5.670374419E-8


############### Set Script Inputs
# Primary Radiator Parameters (Conservative Values)
eps = 0.8 # Limit - 0.8
L = 3.5 # m, Limit - 4 m
W = 1.5*4 # m, Limit - 2 m

# He-Xe Coolant Molar Mass (Conservative Value)
MM = 0.04 # kg / mol, Limit - 0.02 kg / mol

# He-Xe Coolant Ratio of Specific Heats
gam = 5 / 3

# Turbine Pressure Ratio
TPR = 1/2

# Reactor System Pressure Factor
Fp = 1.2

# Typical Generator Efficiency
eta_gen = 0.9

# Number of Calculation Points in Each Input Range
divs = 256

# Reactor Outlet Temperature Range
To = np.linspace(900, 1800, divs) # K

# Reactor Power Range
Qin = np.linspace(100E3, 1000E3, divs) # W


############### Initialise Design Constants and Solution Fields
# Calculate HeXe Cp Value
Cp = 2.5 * (R / MM) # J / (kg * K)

# Calculate Turbine Outlet or Radiator Inlet Temperatures
Tr = (TPR**((gam - 1) / gam)) * To

# Initialise Radiator Outlet Solution Field
Te = np.NaN * np.ones((divs, divs)) # K

# Initialise Reactor Inlet Solution Field
Ti = np.NaN * np.ones((divs, divs)) # K

# Initialise Mass Flow Rate Solution Field
m_dot = np.NaN * np.ones((divs, divs)) # kg / s

# Initialise Thermal Conversion Efficiency Solution Field
eta_th = np.NaN * np.ones((divs, divs))


############### Define Functions
# Define Zero Function to Solve for Reactor Inlet Temperature
def get_T_e(x, q, to, tr):
    return ((((1 / (3 * x**3)) - (1 / (3 * tr**3))) \
             * (1 / (to - (x * (Fp / TPR)**((gam - 1) / gam))))) \
             - ((eps * s_b * L * W) / q))

# Define Helper Function to Solve for Mass Flow Rate
def get_m_dot(q, to, ti):
    return (q / (Cp * (to - ti))) # kg / s

# Define Helper Function to Solve for Thermal Conversion Efficiency
def get_eta(to, tr, te, ti):
    return (eta_gen * ((to - tr - ti + te) / (to - ti)))


############### Calculate Solution Fields
Tol = 1E-12
Ta = 4 # K
To_It = np.nditer(To, flags=['f_index'])
Qin_It = np.nditer(Qin, flags=['f_index'])

for to in To_It:
    tr = Tr[To_It.index]
    for q in Qin_It:
        # Find Radiator Outlet Temperature
        def opt_func(x):
            return (get_T_e(x, q, to, tr))

        try:
            res = root_scalar(opt_func, bracket=[Ta, tr], x0=Ta,
                              xtol=Tol,
                              method='secant')
        except ValueError:
            continue

        # Handle Case for No Convergence by Rejecting the Result
        if res.converged is not True:
            break

        # Save Index Values for Both Iterator to Save Results
        i = To_It.index
        j = Qin_It.index

        te = res.root
        ti = te * ((Fp / TPR)**((gam - 1) / gam))
        Te[i, j] = te
        Ti[i, j] = ti
        m_dot[i, j] = get_m_dot(q, to, ti)
        eta_th[i, j] = get_eta(to, tr, te, ti)

    # Reset Reactor Power Values Iterator for Next Temperature Value
    Qin_It.reset()

    # Clear Console Output and Show Solution Progress
    if name == 'nt':
        system('cls')
    else:
        system('clear')

    print('Solution Fields >> Calculation Progress >> {:2.1f} %'.format(
            ((100 * To_It.index) / To.shape[0])))


############### Plot Solution Fields
###### Reactor Inlet Temperature with Fixed Range for Comparison
fig = plt.figure()
ax = fig.add_subplot(111)
plt_title = 'Reactor Inlet Temperatures at' + \
    ' TPR of {:1.2f}'.format(TPR)
ax.set_title(plt_title)
ax.set_ylabel('Reactor Outlet Temperature (K)')

map = plt.colormaps.get_cmap('plasma')
map.set_bad('black', 1.0)
mesh = ax.pcolormesh(Qin, To, Ti, vmin=400, vmax=1800,
                     cmap=map,
                     shading='nearest',
                     rasterized=True)

# Plot and Configure Colorbar Labels
cb = plt.colorbar(mesh)
cb.ax.set_title('$T_i$ (K)', fontsize=8, loc='left')
cb.ax.tick_params(labelsize=8)

# Update X Axis Values to kW
ax.set_xlabel('Reactor Thermal Power (kW)')
ax.xaxis.set_major_formatter(lambda x, _: str(int(x * 1E-3)))

# Finalise and Save Reactor Inlet Temperature Plot
fig.tight_layout()
fig.savefig(f'{TPR:1.2f}' + '_TRP_Gas_Inlet_Temp.pdf', format='pdf')

# Clear Figure and Axes
fig.clear()
ax.clear()

###### Mass Flow Rate with Fixed Range for Comparison
fig = plt.figure()
ax = fig.add_subplot(111)
plt_title = 'HeXe Mass Flow Rates at' + \
    ' TPR of {:1.2f}'.format(TPR)
ax.set_title(plt_title)
ax.set_xlabel('Reactor Thermal Power (W)')
ax.set_ylabel('Reactor Outlet Temperature (K)')

map = plt.colormaps.get_cmap('plasma')
map.set_bad('black', 1.0)
mesh = ax.pcolormesh(Qin, To, np.log10(m_dot),
                       vmin=-1, vmax=1,
                       cmap=map,
                       shading='nearest',
                       rasterized=True)

# Plot and Configure Colorbar Labels
cb = plt.colorbar(mesh, format=lambda x, _: '{:3.2f}'.format(10**x))
cb.ax.set_title('$\\dot{m}$ (kg/s)', fontsize=8, loc='left')
cb.ax.tick_params(labelsize=8)

# Update X Axis Values to kW
ax.set_xlabel('Reactor Thermal Power (kW)')
ax.xaxis.set_major_formatter(lambda x, _: str(int(x * 1E-3)))

# Finalise and Save Reactor Mass Flow Plot
fig.tight_layout()
fig.savefig(f'{TPR:1.2f}' + '_TRP_Gas_Mass_Flow.pdf', format='pdf')

# Clear Figure and Axes
fig.clear()
ax.clear()

###### System Electrical Conversion Efficiency with Fixed Range for Comparison
fig = plt.figure()
ax = fig.add_subplot(111)
plt_title = 'Electrical Conversion Efficiencies at' + \
    ' TPR of {:1.2f}'.format(TPR)
ax.set_title(plt_title)
ax.set_xlabel('Reactor Thermal Power (W)')
ax.set_ylabel('Reactor Outlet Temperature (K)')

map = plt.colormaps.get_cmap('plasma')
map.set_bad('black', 1.0)
mesh = ax.pcolormesh(Qin, To, eta_th, vmin=0, vmax=0.3,
                       cmap=map,
                       shading='nearest',
                       rasterized=True)

# Plot and Configure Colorbar Labels
cb = plt.colorbar(mesh, format=lambda x, _: '{:.1%}'.format(x))
cb.ax.set_title('$\\eta_{th}$', fontsize=8, loc='left')
cb.ax.tick_params(labelsize=8)

# Update X Axis Values to kW
ax.set_xlabel('Reactor Thermal Power (kW)')
ax.xaxis.set_major_formatter(lambda x, _: str(int(x * 1E-3)))

# Finalise and Save Reactor Conversion Efficiencies Plot
fig.tight_layout()
fig.savefig(f'{TPR:1.2f}' + '_TRP_Gas_Conv_Eta.pdf', format='pdf')

# Clear Figure and Axes
fig.clear()
ax.clear()
