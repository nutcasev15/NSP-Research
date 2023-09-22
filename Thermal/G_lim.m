% Script to Find Bounding Reactor Core Mass Flux and Mass Flow for NSP Systems
% Refer ntrs.nasa.gov./citations./20060056311 for HeXe Coolant Data
% Refer doi.org/10.1016/0009-2614(75)80286-7 for He Coolant Data
% Refer iaea.org/publications/7965/thermophysical-properties-of\
%       -materials-for-nuclear-engineering-a-tutorial-and-collection-of-data
%       for UN Integrated Thermal Conductivity Data
% Refer osti.gov/biblio/532973 for MA956 ODS Steel Thermal Conductivity Data

clc;
clear all;
close all;


% ########################### Define Natural Constants
% Import Physical Constants in SI Units
R = 8.314462618;


% ########################### Set Script Inputs
% Coolant Inlet Temperature (Minimum Possible)
Ti = 400; % K

% Coolant Outlet Temperature (Maximum Possible)
To = 1200; % K

% Fuel Element Linear Power Rating (Maximum Possible)
q_dash = 275E2; % W / m

% Number of Fuel Element Rings in Reactor Core
Rings = 8;

% Core Height
H = 0.46; % m

% MA 956 ODS Steel Clad Radius
rc = 0.84; % cm

% Fuel Element Coolant Channel Area
A = 0.75; % cm^2

% Coolant Nominal Pressure
P = 2E6; % Pa

% Coolant Molar Mass
MM = 40; % g / mol

% Coolant Ratio of Specific Heats
gam = 5 / 3;

% K Value for Pressure Drop Limit
Kf = 1.5;

% Core Pressure Drop Limit
dP = 0.2E6; % Pa


% ########################### Import Coolant Thermal Data
% Import HeXe Coolant Data from Johnson-2006
mu_Data = csvread('HeXe_mu.csv');
mu = @ (T) (interp2(mu_Data(1, 2:end),
                        mu_Data(2:end, 1),...
                        mu_Data(2:end, 2:end),...
                        MM, T));

##% Import He Coolant Data from Jain-1975
##mu_Data = csvread('He_mu.csv');
##mu = @ (T) (interp1(mu_Data(:, 1), mu_Data(:, 2), T));


% ########################### Calculate Hydraulic Diameter of Coolant Channel
% Find Channel Outer Radius
ro = sqrt((A ./ pi) + rc .^ 2); % cm

% Calculate Hydraulic Diameter
D_h = 2 .* (ro - rc) .* 1E-2; % m


% ########################### Calculate Limiting Mass Flux Values
% Calculate Incompressible Flow Mass Flux Limit
G_Mach = (0.3 .* P .* sqrt((gam .* MM .* 1E-3) ./ (R .* To))); % Kg / (m^2 * s)

% Calculate Taylor-1988 HTC Correlation Reynolds Number Mass Flux Limit
G_HTC = (6E4 .* mu(Ti)) ./ D_h; % Kg / (m^2 * s)

% Find Inlet and Outlet Coolant Density at Nominal Pressure
rho_i = P ./ (R ./ (MM .* 1E-3) .* Ti); % Kg / m^3
rho_o = P ./ (R ./ (MM .* 1E-3) .* To); % Kg / m^3

% Calculate Darcy-Weishbach Pressure Loss Mass Flux Limit
% Use Blasius Relation for f Factor
% Assume Maximum Reynolds Number of Taylor-1988 HTC Correlation
G_dP = sqrt(dP ./ (((0.316 .* (6E4 .^ -0.25)) .* H .* Kf ./ (2 .* rho_o .* D_h))...
                  + 0.5 .* (rho_o .^ -1 - rho_i .^ -1))); % Kg / (m^2 * s)

% Calculate Basic Thermal Cooling Mass Flux Limit
G_TH = q_dash .* H...
       ./ ((2.5 .* R ./ (MM .* 1E-3)) .* A .* 1E-4 .* (To - Ti)); % Kg / (m^2 * s)


% ########################### Calculate Corresponding Core Mass Flow Values
% Find Number of Elements in Reactor Core
FE_Num = ((3 .* (Rings .^ 2)) - (3 .* Rings) + 1);

% Calculate Limiting Core Mass Flow for Incompressible Flow
M_Mach = G_Mach .* FE_Num .* A .* 1E-4; % Kg / s

% Calculate Limiting Core Mass Flow for Taylor-1988 HTC Correlation
M_HTC = G_HTC .* FE_Num .* A .* 1E-4; % Kg / s

% Calculate Limiting Core Mass Flow for Maximum Core Pressure Loss Limit
M_dP = G_dP .* FE_Num .* A .* 1E-4; % Kg / s

% Calculate Minimum Core Mass Flow for Input Reactor Core Thermal Power
M_TH = G_TH .* FE_Num .* A .* 1E-4; % Kg / s


% ########################### Output Calculated Values
clc;
printf('\nScript to Calculate Bounding Reactor Core Mass Flux for NSP Systems');
printf('\n\n');
printf('######## Script Inputs: \n');
printf('%-60s %20.3E K\n', 'Coolant Inlet Temperature:', Ti);
printf('%-60s %20.3E K\n', 'Coolant Outlet Temperature:', To);
printf('%-60s %20.3E W/m\n', 'Maximum Fuel Linear Power Rating:', q_dash);
printf('%-60s %20.3E\n', 'Number of Fuel Element Rings:', Rings);
printf('%-60s %20.3E m\n', 'Reactor Core Height:', H);
printf('\n%-60s %20.3E W\n\n', 'Reactor Core Maximum Thermal Power:', FE_Num .* q_dash .* H);
printf('%-60s %20.3E cm\n', 'Fuel Element Clad Radius:', rc);
printf('%-60s %20.3E cm^2\n', 'Fuel Element Coolant Channel Area:', A);
printf('%-60s %20.3E Pa\n', 'Nominal Coolant Pressure:', P);
printf('%-60s %20.3E g/mol\n', 'Coolant Molar Mass:', MM);
printf('%-60s %20.3E\n', 'Coolant Ratio of Specific Heats:', gam);
printf('%-60s %20.3E\n', 'K Value for Core Pressure Drop Limit:', Kf);
printf('%-60s %20.3E Pa\n', 'Core Pressure Drop Limit:', dP);

printf('\n######## Limiting Mass Flux Values\n');
printf('%-60s %20.3E Kg/(m^2*s)\n', 'Incompressible Flow Core Mass Flux Limit:', G_Mach);
printf('%-60s %20.3E Kg/(m^2*s)\n', 'Taylor HTC Correlation Mass Flux Limit:', G_HTC);
printf('%-60s %20.3E Kg/(m^2*s)\n', 'Darcy-Weishbach Pressure Loss Mass Flux Limit:', G_dP);
printf('%-60s %20.3E Kg/(m^2*s)\n', 'Fundamental Core Cooling Minimum Mass Flux Limit:', G_TH);

printf('\n######## Limiting Reactor Core Mass Flow Rates\n');
printf('%-60s %20.3E Kg/s\n', 'Incompressible Flow Core Mass Flow Rate Limit:', M_Mach);
printf('%-60s %20.3E Kg/s\n', 'Taylor HTC Correlation Mass Flow Rate Limit:', M_HTC);
printf('%-60s %20.3E Kg/s\n', 'Darcy-Weishbach Pressure Loss Mass Flow Rate Limit:', M_dP);
printf('%-60s %20.3E Kg/s\n', 'Fundamental Core Cooling Minimum Mass Flow Rate Limit:', M_TH);

