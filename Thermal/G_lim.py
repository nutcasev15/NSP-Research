# Bounding Reactor Core Mass Flux and Mass Flow Calculator for NSP Systems
# Powered by NumPy 1.26.4 and SciPy 1.14.0
#
# Refer ntrs.nasa.gov/citations/20060056311 for HeXe Coolant Data
# Refer doi.org/10.1016/0009-2614(75)80286-7 for He Coolant Data
# Refer iaea.org/publications/7965/thermophysical-properties-of\
#       -materials-for-nuclear-engineering-a-tutorial-and-collection-of-data
#       for UN Integrated Thermal Conductivity Data
# Refer osti.gov/biblio/532973 for MA956 ODS Steel Thermal Conductivity Data


############### Library Imports
import numpy as np
from scipy.interpolate import RegularGridInterpolator, CubicSpline


############### Define Natural Constants
# Import Physical Constants in SI Units
R = 8.314462618


############### Set Script Inputs
# Coolant Inlet Temperature (Minimum Possible)
Ti = 400 # K

# Coolant Outlet Temperature (Maximum Possible)
To = 1200 # K

# Fuel Element Linear Power Rating (Maximum Possible)
Q = 250E2 # W/m

# Number of Fuel Element Rings in Reactor Core
Rings = 5

# Core Height (Maximum Possible)
H = 0.60 # m

# MA956 ODS Steel Clad Radius
Rc = 0.89 # cm

# Fuel Element Coolant Channel Area
A = 0.75 # cm^2

# Coolant Nominal Pressure
P = 2E6 # Pa

# Coolant Molar Mass
MM = 40 # g / mol

# Coolant Ratio of Specific Heats
gam = 5 / 3

# K Value for Pressure Drop Limit
Kf = 1.5

# Core Pressure Drop Limit
dP = 0.2E6 # Pa


############### Import Coolant Thermal Data
# Import HeXe Coolant Data from Johnson-2006
CSVData = np.genfromtxt('HeXe_mu.csv', delimiter=',')
muIP = RegularGridInterpolator((CSVData[1:, 0], CSVData[0, 1:]),
                               CSVData[1:, 1:], 'cubic')

# Import He Coolant Data from Jain-1975
# CSVData = np.genfromtxt('He_mu.csv', delimiter=',')
# muIP = CubicSpline(CSVData[:, 0], CSVData[:, 1])

# Define Function to Return Viscosity Data vs Temperature
def mu(T):
    try:
        # Succeeds with Jain-1975 He Dataset
        return muIP(T)
    except IndexError:
        # Succeeds with Johnson-2006 HeXe Dataset
        return muIP((T, MM))


############### Calculate Hydraulic Diameter of Coolant Channel
# Find Channel Outer Radius
Ro = np.sqrt((A / np.pi) + Rc**2) # cm

# Calculate Hydraulic Diameter
Dh = 2 * (Ro - Rc) * 1E-2 # m


############### Calculate Limiting Mass Flux Values
# Calculate Incompressible Flow Mass Flux Limit
G_Mach = (0.3 * P * np.sqrt((gam * MM * 1E-3) / (R * To))) # kg / (m^2 * s)

# Calculate Taylor-1988 HTC Correlation Reynolds Number Mass Flux Limit
G_HTC = (6E4 * mu(Ti)) / Dh # kg / (m^2 * s)

# Find Inlet and Outlet Coolant Density at Nominal Pressure
Rho_IP = P / (R / (MM * 1E-3) * Ti) # kg / m^3
Rho_OP = P / (R / (MM * 1E-3) * To) # kg / m^3

# Calculate Darcy-Weishbach Pressure Loss Mass Flux Limit
# Use Blasius Relation for f Factor
# Assume Maximum Reynolds Number of Taylor-1988 HTC Correlation
G_dP = np.sqrt(dP / (((0.316 * (6E4**-0.25)) * H * Kf / (2 * Rho_OP * Dh)) \
                     + 0.5 * (Rho_OP**-1 - Rho_IP**-1))) # kg / (m^2 * s)

# Calculate Basic Thermal Cooling Mass Flux Limit
G_TH = Q * H / ((2.5 * R / (MM * 1E-3)) \
                * A * 1E-4 * (To - Ti)) # kg / (m^2 * s)


############### Calculate Corresponding Core Mass Flow Values
# Find Number of Elements in Reactor Core
FE_Num = ((3 * (Rings**2)) - (3 * Rings) + 1)

# Calculate Limiting Core Mass Flow for Incompressible Flow
M_Mach = G_Mach * FE_Num * A * 1E-4 # kg / s

# Calculate Limiting Core Mass Flow for Taylor-1988 HTC Correlation
M_HTC = G_HTC * FE_Num * A * 1E-4 # kg / s

# Calculate Limiting Core Mass Flow for Maximum Core Pressure Loss Limit
M_dP = G_dP * FE_Num * A * 1E-4 # kg / s

# Calculate Minimum Core Mass Flow for Input Reactor Core Thermal Power
M_TH = G_TH * FE_Num * A * 1E-4 # kg / s


############### Output Calculated Values
print(' ')
print('Script to Calculate Bounding Reactor Core Mass Flux for NSP Systems')
print(' ')
print('######## Script Inputs:\n')
print('{:<60} {:20.3f} K'.format('Coolant Inlet Temperature:', Ti))
print('{:60} {:20.3f} K'.format('Coolant Outlet Temperature:', To))
print('{:60} {:20.3f} W/m'.format('Maximum Fuel Linear Power Rating:', Q))
print('{:60} {:20.3f}'.format('Number of Fuel Element Rings:', Rings))
print('{:60} {:20.3f} m'.format('Reactor Core Height:', H))
print('\n{:60} {:20.3f} W\n'.format('Reactor Core Maximum Thermal Power:', \
                                      FE_Num * Q * H))
print('{:60} {:20.3f} cm'.format('Fuel Element Clad Radius:', Rc))
print('{:60} {:20.3f} cm^2'.format('Fuel Element Coolant Channel Area:', A))
print('{:60} {:20.3f} Pa'.format('Nominal Coolant Pressure:', P))
print('{:60} {:20.3f} g/mol'.format('Coolant Molar Mass:', MM))
print('{:60} {:20.3f}'.format('Coolant Ratio of Specific Heats:', gam))
print('{:60} {:20.3f}'.format('K Value for Core Pressure Drop Limit:', Kf))
print('{:60} {:20.3f} Pa\n'.format('Core Pressure Drop Limit:', dP))

print('######## Limiting Mass Flux Values\n')
print('{:60} {:20.3f} kg/m^2*s'.format(
        'Incompressible Flow Core Mass Flux Limit:', G_Mach))
print('{:60} {:20.3f} kg/m^2*s'.format(
        'Taylor HTC Correlation Mass Flux Limit:', G_HTC))
print('{:60} {:20.3f} kg/m^2*s'.format(
        'Darcy-Weishbach Pressure Loss Mass Flux Limit:', G_dP))
print('{:60} {:20.3f} kg/m^2*s\n'.format(
        'Fundamental Core Cooling Minimum Mass Flux Limit:', G_TH))

print('######## Limiting Reactor Core Mass Flow Rates\n')
print('{:60} {:20.3f} kg/s'.format(
        'Incompressible Flow Core Mass Flow Rate Limit:', M_Mach))
print('{:60} {:20.3f} kg/s'.format(
        'Taylor HTC Correlation Mass Flow Rate Limit:', M_HTC))
print('{:60} {:20.3f} kg/s'.format(
        'Darcy-Weishbach Pressure Loss Mass Flow Rate Limit:', M_dP))
print('{:60} {:20.3f} kg/s\n'.format(
        'Fundamental Core Cooling Minimum Mass Flow Rate Limit:', M_TH))
