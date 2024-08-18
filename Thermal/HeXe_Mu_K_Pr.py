# Transport Properties Generator for HeXe Mixtures for NSP Systems
# Powered by NumPy 1.26.4 and SciPy 1.14.0
#
# Refer ntrs.nasa.gov/citations/20060056311 for Equations and Validation
# Refer physics.nist.gov/cgi-bin/Compositions/stand_alone.pl for Molar Masses


############### Library Imports
import numpy as np
from scipy.interpolate import RegularGridInterpolator, CubicSpline


############### Define Natural Constants
# Import Physical Constants in SI Units
R = 8.314462618


############### Define Lennard-Jones Constants
# Sigma, Reduced Temperature and Molar Mass for He
s1 = 2.576 # A
ek1 = 10.22 # K
M1 = 4.002 # g / mol

# Sigma, Reduced Temperature and Molar Mass for Xe
s2 = 4.055 # A
ek2 = 229 # K
M2 = 131.293 # g / mol


############### Set Script Inputs
T = np.arange(400, 2100, 100) # K
MM = np.arange(10, 90, 5) # g / mol

# Duplicate Array Data for Calculation Vectorisation in 2D Arrays
T = np.tile(T, (MM.shape[0], 1))
MM = np.tile(MM, (T.shape[1], 1))
T = T.transpose()


############### Import Data from Johnson-2006
# Collision Integral Evaluation vs Reduced Temperature
CIV_Data = np.genfromtxt('CIV.csv', delimiter=',')
CIV = CubicSpline(CIV_Data[:, 0], CIV_Data[:, 1])

# 1st Lennard-Jones Potential Function vs Reduced Temperature
ALJ_Data = np.genfromtxt('ALJ.csv', delimiter=',')
ALJ = CubicSpline(ALJ_Data[:, 0], ALJ_Data[:, 1])

# 2nd Lennard-Jones Potential Function vs Reduced Temperature
BLJ_Data = np.genfromtxt('BLJ.csv', delimiter=',')
BLJ = CubicSpline(BLJ_Data[:, 0], BLJ_Data[:, 1])

# 3rd Order Correction Coefficients for Thermal Conductivity from Singh-1992
flm_Data = np.genfromtxt('flm.csv', delimiter=',')
flm = RegularGridInterpolator((flm_Data[1:, 0], flm_Data[0, 1:]),
                              flm_Data[1:, 1:], 'cubic')


############### Calculate Initial Calculation Variables
# Mixture Specific Heat
Cp = 2.5 * (R / (MM * 1E-3)) # J / (kg * K)

# Mixture Molar Fractions for He and Xe
x1 = (M2 - MM) / (M2 - M1)
x2 = 1 - x1

# Lennard-Jones Parameters for Equivalent Single Component Gas
s12 = (1 / 2) * (s1 + s2) # A
ek12 = np.sqrt(ek1 * ek2) # K

# Reduced Temperatures for He, Xe and Equivalent Single Component Gas
T1 = T / ek1
T2 = T / ek2
T12 = T / ek12

# Viscosity for He, Xe and Equivalent Single Component Gas
mu1 = 266.93 * 1E-7 * (np.sqrt(M1 * T) / (s1**2 * CIV(T1))) # g / (cm * s)
mu2 = 266.93 * 1E-7 * (np.sqrt(M2 * T) / (s2**2 * CIV(T2))) # g / (cm * s)
mu12 = 266.93 * 1E-7 * (np.sqrt(2 * M1 * M2 * T / (M1 + M2)) \
                         / (s12**2 * CIV(T12))) # g / (cm * s)

# Thermal Conductivity for He, Xe and Equivalent Single Component Gas
k1 = 1989.1 * 1E-7 * (np.sqrt(T / M1) / (s1**2 * CIV(T1))) # cal / (cm * s * K)
k2 = 1989.1 * 1E-7 * (np.sqrt(T / M2) / (s2**2 * CIV(T2))) # cal / (cm * s * K)
k12 = 1989.1 * 1E-7 * (np.sqrt(T * (M1 + M2) / (2 * M1 * M2)) \
                        / (s12**2 * CIV(T12))) # cal / (cm * s * K)


############### Calculate Viscosity of the HeXe Mixtures
# Calculate Intermediate Variables X, Y and Z
Xmu = (x1**2 / mu1) + ((2 * x1 * x2) / mu12) + (x2**2 / mu2)
Ymu = (3 * ALJ(T12) / 5) * ((x1**2 / mu1) * (M1 / M2)  \
      + ((2 * x1 * x2) / mu12) \
        * ((M1 + M2)**2 / (4 * M1 * M2)) * (mu12**2 / (mu1 * mu2)) \
      + (x2**2 / mu2) * (M2 / M1))
Zmu = (3 * ALJ(T12) / 5) * ((x1**2) * (M1 / M2) \
      + (2 * x1 * x2) * (((M1 + M2)**2 / (4 * M1 * M2)) \
        * ((mu12 / mu1) + (mu12 / mu2)) - 1) \
      + (x2**2) * (M2 / M1))

# Calculate Mixture Viscosity
mu = (1 + Zmu) / (Xmu + Ymu)  # g / (cm * s)

# Convert Viscosity to SI Units
mu *= 0.1


############### Calculate Thermal Conductivity of the HeXe Mixtures
# Calculate U Parameters for X, Y and Z
U1 = (4 / 15) * ALJ(T12) \
     - (1 / 12) * ((12 / 5) * BLJ(T12) + 1) * (M1 / M2) \
     + (1 / 2) * ((M1 - M2)**2 / (M1 * M2))
U2 = (4 / 15) * ALJ(T12) \
     - (1 / 12) * ((12 / 5) * BLJ(T12) + 1) * (M2 / M1) \
     + (1 / 2) * ((M2 - M1)**2 / (M1 * M2))
UY = (4 / 15) * ALJ(T12) * ((M1 + M2)**2 / (4 * M1 * M2)) \
               * (k12**2 / (k1 * k2)) \
     - (1 / 12) * ((12 / 5) * BLJ(T12) + 1) \
     - (5 / (32 * ALJ(T12))) * ((12 / 5) * BLJ(T12) - 5) \
                               * ((M1 - M2)**2 / (M1 * M2))
UZ = (4 / 15) * ALJ(T12) * (((M1 + M2)**2 / (4 * M1 * M2)) \
               * ((k12 / k1) + (k12 / k2)) - 1) \
     - (1 / 12) * ((12 / 5) * BLJ(T12) + 1)

# Calculate Intermediate Variables X, Y and Z using Calculated U Parameters
Xlm = (x1**2 / k1) + ((2 * x1 * x2) / k12) + (x2**2 / k2)
Ylm = (x1**2 / k1) * U1 + ((2 * x1 * x2) / k12) * UY + (x2**2 / k2) * U2
Zlm = (x1**2) * U1 + (2 * x1 * x2) * UZ + (x2**2) * U2

# Calculate Mixture Thermal Conductivity
k3 = (1 + Zlm) / (Xlm + Ylm) # cal / (cm * s * K)

# Convert Thermal Conductivity to SI Units
k3 *= 418.4

# Apply Singh's 3rd Order Thermal Conductivity Correction
k = k3 * flm((x1, T))


############### Calculate Prandtl Number of the HeXe Mixtures
Pr = Cp * mu / k


############### Assemble and Save CSV Files
# Create CSV Data Container
CSVContainer = np.zeros((T.shape[0] + 1, MM.shape[1] + 1))
CSVContainer[1:, 0] = T[:, 0]
CSVContainer[0, 1:] = MM[0, :]

# Load Negative Placeholder Value in 1st Cell
CSVContainer[0, 0] = -1

# Save Viscosity Data to CSV
CSVContainer[1:, 1:] = mu # kg / (m * s)
np.savetxt('HeXe_mu.csv', CSVContainer, delimiter=',', fmt='%G')

# Save Thermal Conductivity Data to CSV
CSVContainer[1:, 1:] = k # W / (m * K)
np.savetxt('HeXe_k.csv', CSVContainer, delimiter=',', fmt='%G')

# Save Prandtl Number Data to CSV
CSVContainer[1:, 1:] = Pr
np.savetxt('HeXe_Pr.csv', CSVContainer, delimiter=',', fmt='%G')
