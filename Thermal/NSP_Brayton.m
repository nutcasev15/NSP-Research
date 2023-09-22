% Script to Find Bounding Thermodynamic Parameters for Simple NSP Systems
% Refer NSP_Turbine.drawio for Equations
% Refer ec.europa.eu/research/participants/documents/downloadPublic?\
%       documentIds=080166e5c802fe77&appId=PPGMS for Radiator Parameters
% Refer doi.org/10.2514/1.27664 for He-Xe Coolant Data

clc;
clear all;
close all;


% ########################### Define Natural Constants
% Import Physical Constants in SI Units
R = 8.314462618;
s_b = 5.670374419E-8;


% ########################### Set Script Inputs
% Primary Radiator Parameters (Conservative Values)
epsilon = 0.8; % Limit - 0.8
L = 3.5; % m, Limit - 4 m
W = 1.5*4; % m, Limit - 2 m

% He-Xe Coolant Molar Mass (Conservative Value)
MM = 0.04; % kg / mol, Limit - 0.02 kg / mol

% He-Xe Coolant Ratio of Specific Heats
gam = 5 / 3;

% Turbine Pressure Ratio
TPR = 1/2;

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


% ########################### Initialise Design Constants and Solution Fields
% Calculate HeXe Cp Value
C_p = 2.5 * (R / MM); % J / (kg * K)

% Calculate Turbine Outlet or Radiator Inlet Temperatures
T_r = (TPR.^((gam - 1) ./ gam)) .* T_o;

% Initialise Radiator Outlet Solution Field
T_e = NaN(divs, divs); % K

% Initialise Reactor Inlet Solution Field
T_i = NaN(divs, divs); % K

% Initialise Mass Flow Rate Solution Field
m_dot = NaN(divs, divs); % kg / s

% Initialise Thermal Conversion Efficiency Solution Field
eta_th = NaN(divs, divs);


% ########################### Define Lambda Functions
% Define Zero Function to Solve for Reactor Inlet Temperature
T_e_func = @ (x, Q_in, T_o, T_r) ((((1 ./ (3 * x.^3)) - (1 ./ (3 .* T_r.^3))) ...
           .* (1 ./ (T_o - (x .* (F_p ./ TPR).^((gam - 1) ./ gam)))))...
           - ((epsilon .* s_b .* L .* W) ./ Q_in));

% Define Helper Function to Solve for Mass Flow Rate
m_func = @ (Q_in, T_o, T_i) (Q_in ./ (C_p .* (T_o - T_i))); % kg / s

% Define Helper Function to Solve for Thermal Conversion Efficiency
eta_func = @ (T_o, T_r, T_e, T_i) (eta_gen .* ((T_o - T_r - T_i + T_e)...
                                              ./ (T_o - T_i)));


% ########################### Calculate Solution Fields
Tol = 1E-12;
T_a = 4; % K
for i = 1:divs
  t_o = T_o(i);
  t_r = T_r(i);
  for j = 1:divs
    Q = Q_in(j);

    try
      % Find Radiator Outlet Temperature
      tmp_func = @ (x) (T_e_func(x, Q, t_o, t_r));
      [t_e, val, info, out] = fzero(tmp_func, [T_a t_r],...
                                    optimset('TolX', Tol, 'Display', 0));
    catch
      break;
    end_try_catch

    % Handle Case for Singular Point Convergence by Rejecting the t_e Value
    if (info == -5)
      continue;
    endif

    t_i = t_e .* ((F_p ./ TPR).^((gam - 1) ./ gam));
    T_e(i, j) = t_e;
    T_i(i, j) = t_i;
    m_dot(i, j) = m_func(Q, t_o, t_i);
    eta_th(i, j) = eta_func(t_o, t_r, t_e, t_i);
  endfor
  clc;
  printf('Solution Fields >> Calculation Progress >> %2.1f %%\n', ((100 * i) / divs));
endfor
clc;


% ########################### Plot Solution Fields
% Plot for Reactor Inlet Temperature with Fixed Range for Comparison
figure(1, 'name', 'Plot of Reactor Inlet Temperature');
imagesc(Q_in, T_o, T_i, [400 1800]);
set(gca, 'ydir', 'normal');
colormap(cubehelix(divs));

cb = colorbar();
labels = linspace(400, 1800, 32);
set(cb, 'ytick', labels);
set(cb, 'yticklabel', strsplit(num2str(labels, '%.0f ')));
set(cb, 'title', 'T_i (K)');
set(cb, 'fontsize', 10);

str = sprintf('\t\tReactor Inlet Temperatures for Turbine Pressure Ratio of %1.2f\n', TPR);
##title(str);
xlabel('Reactor Thermal Power (W)');
ylabel('Reactor Outlet Temperature (K)');
set(gca, 'fontsize', 12);

% Plot for Mass Flow Rate with Fixed Range for Comparison
figure(2, 'name', 'Logarithmic Plot of Mass Flow Rate');
imagesc(Q_in, T_o, log10(m_dot), [-1 1]);
set(gca, 'ydir', 'normal');
colormap(cubehelix(divs));

cb = colorbar();
labels = linspace(-1, 1, 32);
set(cb, 'ytick', labels);
set(cb, 'yticklabel', strsplit(num2str(power(10, labels), '%0.2f ')));
set(cb, 'title', 'ṁ (Kg/s)');
set(cb, 'fontsize', 10);

str = sprintf('\t\tMass Flow Rates for Turbine Pressure Ratio of %1.2f\n', TPR);
##title(str);
xlabel('Reactor Thermal Power (W)');
ylabel('Reactor Outlet Temperature (K)');
set(gca, 'fontsize', 12);

% Plot for System Electrical Conversion Efficiency with Fixed Range for Comparison
figure(3, 'name', 'System Electrical Conversion Efficiency');
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
clear -x Q_in T_o T_r T_e T_i m_dot eta_th...
         C_p epsilon L W MM gam TPR F_p eta_gen;

