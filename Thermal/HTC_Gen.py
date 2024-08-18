# Heat Transfer Coefficient Table Generator for NSP Systems
# Powered by NumPy 1.26.4 and SciPy 1.14.0
#
# Refer ntrs.nasa.gov/citations/20060056311 for HeXe Coolant Data
# Refer doi.org/10.1016/0009-2614(75)80286-7 for He Coolant Data
# Refer physics.nist.gov/cgi-bin/Compositions/stand_alone.pl for Molar Masses


############### Library Imports
import numpy as np
from scipy.interpolate import RegularGridInterpolator, CubicSpline


############### Set Script Inputs
G = np.arange(100, 1600, 100) # kg / (m^2 * s)
T = np.arange(400, 2100, 100) # K
HeXe_MM = 40 # g / mol
Dh = 0.26354E-2 # m

# Duplicate Array Data for Calculation Vectorisation in 2D Arrays
G = np.tile(G, (T.shape[0], 1))
T  = np.tile(T, (G.shape[1], 1))
T = T.transpose()
HeXe_MM = HeXe_MM * np.ones(G.shape)
Dh = Dh * np.ones(T.shape)


############### Import Processed HeXe Data from Johnson-2006
# Thermal Conductivity vs Temperature at Input Molar Mass
HeXe_k_Data = np.genfromtxt('HeXe_k.csv', delimiter=',')
get_HeXe_k = RegularGridInterpolator((HeXe_k_Data[1:, 0],
                                      HeXe_k_Data[0, 1:]),
                                      HeXe_k_Data[1:, 1:],
                                      'cubic')

# Viscosity vs Temperature at Input Molar Mass
HeXe_mu_Data = np.genfromtxt('HeXe_mu.csv', delimiter=',')
get_HeXe_mu = RegularGridInterpolator((HeXe_mu_Data[1:, 0],
                                       HeXe_mu_Data[0, 1:]),
                                       HeXe_mu_Data[1:, 1:],
                                       'cubic')

# Prandtl Number vs Temperature at Input Molar Mass
HeXe_Pr_Data = np.genfromtxt('HeXe_Pr.csv', delimiter=',')
get_HeXe_Pr = RegularGridInterpolator((HeXe_Pr_Data[1:, 0],
                                       HeXe_Pr_Data[0, 1:]),
                                       HeXe_Pr_Data[1:, 1:],
                                       'cubic')


############### Import He Data from Jain-1975
# Thermal Conductivity vs Temperature
He_k_Data = np.genfromtxt('He_k.csv', delimiter=',')
get_He_k = CubicSpline(He_k_Data[:, 0], He_k_Data[:, 1])

# Viscosity vs Temperature
He_mu_Data = np.genfromtxt('He_mu.csv', delimiter=',')
get_He_mu = CubicSpline(He_mu_Data[:, 0], He_mu_Data[:, 1])

# Prandtl Number vs Temperature
He_Pr_Data = np.genfromtxt('He_Pr.csv', delimiter=',')
get_He_Pr = CubicSpline(He_Pr_Data[:, 0], He_Pr_Data[:, 1])


############### Calculate HTC Using Taylor-1988 for HeXe Mixtures
# Calculate Flow Reynolds Number
HeXe_Re = (G * Dh / get_HeXe_mu((T, HeXe_MM)))

# Enforce Reynolds Number Limits of HTC Correlation in Taylor-1988
HeXe_Re[(HeXe_Re < 1.8E4) | (HeXe_Re > 6E4)] = np.NaN

# Calculate HTC Values with T_W/T_b equal to 1
HeXe_HTC = (get_HeXe_k((T, HeXe_MM)) / Dh) * 0.023 * (HeXe_Re**0.8) \
            * (get_HeXe_Pr((T, HeXe_MM))**0.65)


############### Calculate HTC Using Taylor-1988 for He
# Calculate Flow Reynolds Number
He_Re = (G * Dh / get_He_mu(T))

# Enforce Reynolds Number Limits of HTC Correlation in Taylor-1988
He_Re[(He_Re < 1.8E4) | (He_Re > 6E4)] = np.NaN

# Calculate HTC Values with T_W/T_b equal to 1
He_HTC = (get_He_k(T) / Dh) * 0.023 * (He_Re**0.8) * (get_He_Pr(T)**0.65)


############### Assemble and Save CSV Files
# Create CSV Data Container
CSVContainer = np.zeros((T.shape[0] + 1, G.shape[1] + 1))
CSVContainer[1:, 0] = T[:, 0]
CSVContainer[0, 1:] = G[0, :]

# Save Hydraulic Diameter Input for Reference
# Negate Value for More Visibility in the CSV Files
CSVContainer[0, 0] = -Dh[0, 0] # m

# Save HeXe HTC Data to CSV
CSVContainer[1:, 1:] = HeXe_HTC # W / m^2
np.savetxt('HeXe_HTC.csv', CSVContainer, delimiter=',', fmt='%G')

# Save He HTC Data to CSV
CSVContainer[1:, 1:] = He_HTC # W / m^2
np.savetxt('He_HTC.csv', CSVContainer, delimiter=',', fmt='%G')
