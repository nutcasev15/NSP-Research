# Bounding Thermodynamic Parameter Calculator for LHP Cooled NSP Systems
# Powered by NumPy 1.26.4 and SciPy 1.14.0
#
# Refer NSP_Rankine.pdf for Equations
# Refer ec.europa.eu/research/participants/documents/downloadPublic?\
#       documentIds=080166e5c802fe77&appId=PPGMS for Radiator Parameters
# Refer inis.iaea.org/search/search.aspx?orig_q=RN:27004412 for Sodium Data
# Refer osti.gov/biblio/6856038 for Sodium Entropy and Gamma Data


############### Library Imports
from os import system, name
import numpy as np

import matplotlib.pyplot as plt


############### Define Natural Constants
# Import Physical Constants in SI Units
s_b = 5.670374419E-8


############### Set Script Inputs
# Primary Radiator Reference Parameters (Conservative Values)
eps = 0.8 # Limit - 0.8
L = 3.5 # m, Limit - 4 m
W = 1.5*4 # m, Limit - 2 m

# Turbine Pressure Ratio
TPR = 1/2

# Average Ratio of Specific Heats (900 < T < 1800)
gam = 1.5757

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


############### Initialise Solution Fields
# Initialise Turbine Inlet Pressure Solution Field
Po = np.NaN * np.ones((divs, divs)) # Pa

# Initialise Radiator Inlet Temperature Solution Field
Tr = np.NaN * np.ones((divs, divs)) # K

# Initialise Radiator Inlet Vapour Quality Solution Field
Xr = np.NaN * np.ones((divs, divs))

# Initialise Mass Flow Rate Solution Field
m_dot = np.NaN * np.ones((divs, divs)) # kg / s

# Initialise Normalised Primary Radiator Area Solution Field
A_norm = np.NaN * np.ones((divs, divs))

# Initialise Thermal Conversion Efficiency Solution Field
eta_th = np.NaN * np.ones((divs, divs))


############### Import Sodium Thermophysical Property Data
# Define Liquid Sodium Enthalpy Function (380 < T < 2000)
def get_H_l(T):
    return ((-365.77 \
             + 1.6582 * T \
             - 4.2395E-4 * (T**2) \
             + 1.4847E-7 * (T**3) \
             + 2992.6 / T) * 1E3) # J / kg

# Define Sodium Enthalpy of Vaporisation Function (380 < T < 2000)
def get_H_lg(T):
    return ((393.37 * (1 - (T / 2503.7)) \
             + 4393.6 * ((1 - (T / 2503.7))**0.29302)) * 1E3) # J / kg

# Define Gaseous Sodium Enthalpy Function (380 < T < 2000)
def get_H_g(T):
    return (get_H_l(T) + get_H_lg(T)) # J / kg

# Define Liquid Sodium Entropy Function (380 < T < 2000)
def get_S_l(T):
    return ((7.4402E-1 \
             + 8.2785E-3 * T \
             - 9.1681E-6 * (T**2) \
             + 6.0062E-9 * (T**3) \
             - 2.0202E-12 * (T**4) \
             + 2.7196E-16 * (T**5)) * 1E3) # J / (kg * K)

# Define Gaseous Sodium Entropy Function (380 < T < 2000)
def get_S_g(T):
    return ((2.9487E1 \
             - 6.0974E-2 * T \
             + 7.1262E-5 * (T**2) \
             - 4.3225E-8 * (T**3) \
             + 1.3147E-11 * (T**4) \
             - 1.5835E-15 * (T**5)) * 1E3) # J / (kg * K)

# Define Sodium Saturation Pressure Function (900 < T < 2000)
def get_P_sat(T):
    return (np.exp(11.9463 - 12633.73 / T - 0.4672 * np.log(T)) * 1E6) # Pa

# Define Liquid Sodium Density Function (380 < T < 2000)
def get_Rho_l(T):
    return (219.00 \
            + 275.32 * (1 - T / 2503.7) \
            + 511.58 * ((1 - T / 2503.7)**0.5)) # kg / m^3


############### Define Lambda Functions
# Define Helper Function to Solve for Mass Flow Rate
def get_m_dot(q, to, tr):
    return (q / (get_H_l(to) - get_H_l(tr) + get_H_lg(to))) # kg / s

# Define Helper Function to Solve for Normalised Primary Radiator Area
def get_A_norm(m, xr, tr):
    return ((m * get_H_lg(tr) * xr) \
             / (eps * s_b * (tr**4) * L * W))

# Define Helper Function to Solve for Thermal Conversion Efficiency
def get_eta_th(po, xr, to, tr):
    return (eta_gen \
            * ((get_H_g(to) - (xr * get_H_lg(tr)) - get_H_l(tr) \
            - ((Fp * po * (1 - TPR)) / get_Rho_l(tr))) \
            / (get_H_l(to) - get_H_l(tr) + get_H_lg(to))))


############### Calculate Solution Fields
Tol = 1E-12
To_It = np.nditer(To, flags=['f_index'])
Qin_It = np.nditer(Qin, flags=['f_index'])

for to in To_It:
    tr = Tr[To_It.index]
    for q in Qin_It:
        # Find Turbine Inlet and Outlet Pressures
        po = get_P_sat(to)
        pr = po * TPR

        # Calculate Turbine Outlet Temperature Assuming Constant Gamma
        tr = to * (TPR**((gam - 1) / gam))

        # Calculate Turbine Outlet Vapour Quality Using Entropy Data
        xr = ((get_S_g(to) - get_S_l(tr)) / (get_S_g(tr) - get_S_l(tr)))

        # Save Index Values for Both Iterator to Save Results
        i = To_It.index
        j = Qin_It.index

        Po[i, j] = po
        Tr[i, j] = tr
        Xr[i, j] = xr
        m = get_m_dot(q, to, tr)
        A_norm[i, j] = get_A_norm(m, xr, tr)
        m_dot[i, j] = m
        eta_th[i, j] = get_eta_th(po, xr, to, tr)

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
###### Turbine Inlet Pressure with Fixed Range for Comparison
fig = plt.figure()
ax = fig.add_subplot(111)
plt_title = 'Turbine Inlet Pressures at' + \
    ' TPR of {:1.2f}'.format(TPR)
ax.set_title(plt_title)
ax.set_xlabel('Reactor Thermal Power (W)')
ax.set_ylabel('Reactor Outlet Temperature (K)')

map = plt.colormaps.get_cmap('viridis')
map.set_bad('black', 1.0)
mesh = ax.pcolormesh(Qin, To, np.log10(Po),
                       vmin=4, vmax=7,
                       cmap=map,
                       shading='nearest',
                       rasterized=True)

# Plot and Configure Colorbar Labels
cb = plt.colorbar(mesh, format=lambda x, _: '{:3.2E}'.format(10**x))
cb.ax.set_title('$P_o$ (Pa)', fontsize=8, loc='left')
cb.ax.tick_params(labelsize=8)

# Update X Axis Values to kW
ax.set_xlabel('Reactor Thermal Power (kW)')
ax.xaxis.set_major_formatter(lambda x, _: str(int(x * 1E-3)))

# Finalise and Save Reactor Inlet Pressures Plot
fig.tight_layout()
fig.savefig(f'{TPR:1.2f}' + '_TRP_LHP_Inlet_Pressure.pdf', format='pdf')

# Clear Figure and Axes
fig.clear()
ax.clear()

###### Radiator Inlet Temperature with Fixed Range for Comparison
fig = plt.figure()
ax = fig.add_subplot(111)
plt_title = 'Radiator Inlet Temperatures at' + \
    ' TPR of {:1.2f}'.format(TPR)
ax.set_title(plt_title)
ax.set_ylabel('Reactor Outlet Temperature (K)')

map = plt.colormaps.get_cmap('viridis')
map.set_bad('black', 1.0)
mesh = ax.pcolormesh(Qin, To, Tr, vmin=400, vmax=1800,
                       cmap=map,
                       shading='nearest',
                       rasterized=True)

# Plot and Configure Colorbar Labels
cb = plt.colorbar(mesh)
cb.ax.set_title('$T_r$ (K)', fontsize=8, loc='left')
cb.ax.tick_params(labelsize=8)

# Update X Axis Values to kW
ax.set_xlabel('Reactor Thermal Power (kW)')
ax.xaxis.set_major_formatter(lambda x, _: str(int(x * 1E-3)))

# Finalise and Save Radiator Inlet Temperature Plot
fig.tight_layout()
fig.savefig(f'{TPR:1.2f}' + '_TRP_LHP_Rad_Temp.pdf', format='pdf')

# Clear Figure and Axes
fig.clear()
ax.clear()

###### Mass Flow Rate with Fixed Range for Comparison
fig = plt.figure()
ax = fig.add_subplot(111)
plt_title = 'Sodium Mass Flow Rates at' + \
    ' TPR of {:1.2f}'.format(TPR)
ax.set_title(plt_title)
ax.set_xlabel('Reactor Thermal Power (W)')
ax.set_ylabel('Reactor Outlet Temperature (K)')

map = plt.colormaps.get_cmap('viridis')
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
fig.savefig(f'{TPR:1.2f}' + '_TRP_LHP_Mass_Flow.pdf', format='pdf')

# Clear Figure and Axes
fig.clear()
ax.clear()

###### Normalised Primary Radiator Area with Fixed Range for Comparison
fig = plt.figure()
ax = fig.add_subplot(111)
plt_title = 'Radiator Area vs HeXe System at' + \
    ' TPR of {:1.2f}'.format(TPR)
ax.set_title(plt_title)
ax.set_ylabel('Reactor Outlet Temperature (K)')

map = plt.colormaps.get_cmap('viridis')
map.set_bad('black', 1.0)
mesh = ax.pcolormesh(Qin, To, A_norm, vmin=0, vmax=1.0,
                       cmap=map,
                       shading='nearest',
                       rasterized=True)

# Plot and Configure Colorbar Labels
cb = plt.colorbar(mesh)
cb.ax.set_title('$A_{norm}$', fontsize=8, loc='left')
cb.ax.tick_params(labelsize=8)

# Update X Axis Values to kW
ax.set_xlabel('Reactor Thermal Power (kW)')
ax.xaxis.set_major_formatter(lambda x, _: str(int(x * 1E-3)))

# Finalise and Save Radiator Inlet Vapour Quality Plot
fig.tight_layout()
fig.savefig(f'{TPR:1.2f}' + '_TRP_LHP_A_Norm.pdf', format='pdf')

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

map = plt.colormaps.get_cmap('viridis')
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
fig.savefig(f'{TPR:1.2f}' + '_TRP_LHP_Conv_Eta.pdf', format='pdf')

# Clear Figure and Axes
fig.clear()
ax.clear()
