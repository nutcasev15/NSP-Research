% Script to Find Bounding Thermodynamic Parameters for LHP Cooled NSP Systems
% Refer NSP_LHP.drawio for Equations
% Refer ec.europa.eu/research/participants/documents/downloadPublic?\
%       documentIds=080166e5c802fe77&appId=PPGMS for Radiator Parameters
% Refer inis.iaea.org/search/search.aspx?orig_q=RN:27004412 for Sodium Data
% Refer osti.gov/biblio/6856038 for Sodium Entropy and Gamma Data

clc;
clear all;
close all;


% ########################### Define Natural Constants
% Import Physical Constants in SI Units
s_b = 5.670374419E-8;


% ########################### Set Script Inputs
% Primary Radiator Reference Parameters (Conservative Values)
epsilon = 0.8; % Limit - 0.8
L = 3.5; % m, Limit - 4 m
W = 1.5*4; % m, Limit - 2 m

% Turbine Pressure Ratio
TPR = 1/2;

% Average Ratio of Specific Heats (900 < T < 1800)
gam = 1.5757;

% Reactor System Pressure Factor
F_p = 1.2;

% Typical Generator Efficiency
eta_gen = 0.9;

% Number of Calculation Points in Each Input Range
divs = 256;

% Reactor Outlet Temperature Range
T_o = linspace(900, 1800, divs); % K

% Reactor Power Range
Q_in = linspace(100E3, 1000E3, divs); % W


% ########################### Initialise Solution Fields
% Initialise Turbine Inlet Pressure Solution Field
P_o = NaN(divs, divs); % Pa

% Initialise Radiator Inlet Temperature Solution Field
T_r = NaN(divs, divs); % K

% Initialise Radiator Inlet Vapour Quality Solution Field
X_r = NaN(divs, divs);

% Initialise Mass Flow Rate Solution Field
m_dot = NaN(divs, divs); % Kg / s

% Initialise Normalised Primary Radiator Area Solution Field
A_norm = NaN(divs, divs);

% Initialise Thermal Conversion Efficiency Solution Field
eta_th = NaN(divs, divs);


% ########################### Import Sodium Thermophysical Property Data
% Define Liquid Sodium Enthalpy Function (380 < T < 2000)
H_l = @ (T) ((-365.77...
             + 1.6582 .* T...
             - 4.2395E-4 .* (T.^2)...
             + 1.4847E-7 .* (T.^3)...
             + 2992.6 ./ T) .* 1E3); % J / Kg

% Define Sodium Enthalpy of Vapourisation Function (380 < T < 2000)
H_lg = @ (T) ((393.37 .* (1 - (T ./ 2503.7))...
              + 4393.6 .* ((1 - (T ./ 2503.7)).^0.29302)) .* 1E3); % J / Kg

% Define Gaseous Sodium Enthalpy Function (380 < T < 2000)
H_g = @ (T) (H_l(T) + H_lg(T)); % J / Kg

% Define Liquid Sodium Entropy Function (380 < T < 2000)
S_l = @ (T) ((7.4402E-1...
             + 8.2785E-3 .* T...
             - 9.1681E-6 .* (T.^2)...
             + 6.0062E-9 .* (T.^3)...
             - 2.0202E-12 .* (T.^4)...
             + 2.7196E-16 .* (T.^5)) .* 1E3); % J / (Kg * K)

% Define Gaseous Sodium Entropy Function (380 < T < 2000)
S_g = @ (T) ((2.9487E1...
             - 6.0974E-2 .* T...
             + 7.1262E-5 .* (T.^2)...
             - 4.3225E-8 .* (T.^3)...
             + 1.3147E-11 .* (T.^4)...
             - 1.5835E-15 .* (T.^5)) .* 1E3); % J / (Kg * K)

% Define Sodium Saturation Pressure Function (900 < T < 2000)
P_sat = @ (T) (exp(11.9463 - 12633.73 ./ T - 0.4672 .* log(T)) .* 1E6); % Pa

% Define Liquid Sodium Density Function (380 < T < 2000)
Rho_l = @ (T) (219.00...
              + 275.32 .* (1 - T ./ 2503.7)...
              + 511.58 .* ((1 - T ./ 2503.7).^0.5)); % Kg / m^3


% ########################### Define Lambda Functions
% Define Helper Function to Solve for Mass Flow Rate
m_func = @ (Q_in, T_o, T_r) (Q_in ./ (H_l(T_o)...
                                     - H_l(T_r)...
                                     + H_lg(T_o))); % Kg / s

% Define Helper Function to Solve for Normalised Primary Radiator Area
A_norm_func = @ (m, X_r, T_r) ((m .* H_lg(T_r) .* X_r)...
                         ./ (epsilon .* s_b .* (T_r.^4) .* L .* W));

% Define Helper Function to Solve for Thermal Conversion Efficiency
eta_func = @ (P_o, X_r, T_o, T_r) (eta_gen...
                             .* ((H_g(T_o) - (X_r .* H_lg(T_r)) - H_l(T_r)...
                             - ((F_p .* P_o .* (1 - TPR)) ./ Rho_l(T_r)))...
                             ./ (H_l(T_o) - H_l(T_r) + H_lg(T_o))));


% ########################### Calculate Solution Fields
Tol = 1E-12;
for i = 1:divs
  t_o = T_o(i);
  for j = 1:divs
    Q = Q_in(j);

    % Find Turbine Inlet and Outlet Pressures
    p_o = P_sat(t_o);
    p_r = p_o .* TPR;

    % Calculate Turbine Outlet Temperature Assuming Constant Gamma
    t_r = t_o .* (TPR.^((gam - 1) ./ gam));

    % Calculate Turbine Outlet Vapour Quality Using Entropy Data
    x_r = ((S_g(t_o) - S_l(t_r)) ./ (S_g(t_r) - S_l(t_r)));

    P_o(i, j) = p_o;
    T_r(i, j) = t_r;
    X_r(i, j) = x_r;
    m = m_func(Q, t_o, t_r);
    A_norm(i, j) = A_norm_func(m, x_r, t_r);
    m_dot(i, j) = m;
    eta_th(i, j) = eta_func(p_o, x_r, t_o, t_r);
  endfor
  clc;
  printf('Solution Fields >> Calculation Progress >> %2.1f %%\n', ((100 * i) / divs));
endfor
clc;

% ########################### Plot Solution Fields
% Plot for Turbine Inlet Pressure with Fixed Range for Comparison
figure(1, 'name', 'Plot of Turbine Inlet Pressure');
imagesc(Q_in, T_o, log10(P_o), [4 7]);
set(gca, 'ydir', 'normal');
colormap(cubehelix(divs));

cb = colorbar();
labels = linspace(4, 7, 32);
set(cb, 'ytick', labels);
set(cb, 'yticklabel', strsplit(num2str(power(10, labels), '%.2e ')));
set(cb, 'title', 'P_o (Pa)');
set(cb, 'fontsize', 10);

str = sprintf('\t\tTurbine Inlet Pressures for Turbine Pressure Ratio of %1.2f\n', TPR);
##title(str);
xlabel('Reactor Thermal Power (W)');
ylabel('Reactor Outlet Temperature (K)');
set(gca, 'fontsize', 12);

% Plot for Radiator Inlet Temperature with Fixed Range for Comparison
figure(2, 'name', 'Plot of Radiator Inlet Temperature');
imagesc(Q_in, T_o, T_r, [400 1800]);
set(gca, 'ydir', 'normal');
colormap(cubehelix(divs));

cb = colorbar();
labels = linspace(400, 1800, 32);
set(cb, 'ytick', labels);
set(cb, 'yticklabel', strsplit(num2str(labels, '%.0f ')));
set(cb, 'title', 'T_r (K)');
set(cb, 'fontsize', 10);

str = sprintf('\t\tRadiator Inlet Temperatures for Turbine Pressure Ratio of %1.2f\n', TPR);
##title(str);
xlabel('Reactor Thermal Power (W)');
ylabel('Reactor Outlet Temperature (K)');
set(gca, 'fontsize', 12);

% Plot for Radiator Inlet Vapour Quality with Fixed Range for Comparison
figure(3, 'name', 'Plot of Radiator Inlet Vapour Quality');
imagesc(Q_in, T_o, X_r, [0 1.0]);
set(gca, 'ydir', 'normal');
colormap(cubehelix(divs));

cb = colorbar();
labels = linspace(0, 1.0, 32);
set(cb, 'ytick', labels);
set(cb, 'yticklabel', strsplit(num2str(labels, '%.3f ')));
set(cb, 'title', 'X_r');
set(cb, 'fontsize', 10);

str = sprintf('\t\tRadiator Inlet Vapour Quality for Turbine Pressure Ratio of %1.2f\n', TPR);
##title(str);
xlabel('Reactor Thermal Power (W)');
ylabel('Reactor Outlet Temperature (K)');
set(gca, 'fontsize', 12);

% Plot for Mass Flow Rate with Fixed Range for Comparison
figure(4, 'name', 'Plot of Mass Flow Rate');
imagesc(Q_in, T_o, log10(m_dot), [-1 1]);
set(gca, 'ydir', 'normal');
colormap(cubehelix(divs));

cb = colorbar();
labels = linspace(-1, 1, 32);
set(cb, 'ytick', labels);
set(cb, 'yticklabel', strsplit(num2str(power(10, labels), '%0.3f ')));
set(cb, 'title', 'ṁ (Kg/s)');
set(cb, 'fontsize', 10);

str = sprintf('\t\tMass Flow Rates for Turbine Pressure Ratio of %1.2f\n', TPR);
##title(str);
xlabel('Reactor Thermal Power (W)');
ylabel('Reactor Outlet Temperature (K)');
set(gca, 'fontsize', 12);

% Plot for Normalised Primary Radiator Area with Fixed Range for Comparison
figure(5, 'name', 'Plot of Normalised Primary Radiator Area');
imagesc(Q_in, T_o, A_norm, [0 1.0]);
set(gca, 'ydir', 'normal');
colormap(cubehelix(divs));

cb = colorbar();
labels = linspace(0, 1.0, 32);
set(cb, 'ytick', labels);
set(cb, 'yticklabel', strsplit(num2str(labels, '%0.3f ')));
set(cb, 'title', 'A_{norm}');
set(cb, 'fontsize', 10);

str = sprintf('\t\tNormalised Primary Radiator Area for Turbine Pressure Ratio of %1.2f\n', TPR);
##title(str);
xlabel('Reactor Thermal Power (W)');
ylabel('Reactor Outlet Temperature (K)');
set(gca, 'fontsize', 12);

% Plot for System Electrical Conversion Efficiency with Fixed Range for Comparison
figure(6, 'name', 'Plot of System Electrical Conversion Efficiency');
imagesc(Q_in, T_o, eta_th, [0 0.3]);
set(gca, 'ydir', 'normal');
colormap(cubehelix(divs));

cb = colorbar();
labels = linspace(0, 0.3, 32);
set(cb, 'ytick', labels);
set(cb, 'yticklabel', strsplit(num2str(labels, '%0.3f ')));
set(cb, 'title', '\eta_{th}');
set(cb, 'fontsize', 10);

str = sprintf('\t\tSystem Electrical Conversion Efficiency for Turbine Pressure Ratio of %1.2f\n', TPR);
##title(str);
xlabel('Reactor Thermal Power (W)');
ylabel('Reactor Outlet Temperature (K)');
set(gca, 'fontsize', 12);


% ########################### Clear Unnecessary Variables
clear -x Q_in T_o P_o T_r X_r m_dot A_norm eta_th...
         L W TPR F_p eta_gen;

