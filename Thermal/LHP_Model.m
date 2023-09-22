% Script to Find Bounding Fuel Element Linear Power for Sodium Cooled NSP Systems
% Refer iaea.org/publications/7965/thermophysical-properties-of\
%       -materials-for-nuclear-engineering-a-tutorial-and-collection-of-data
%       for UN Integrated Thermal Conductivity Data
% Refer osti.gov/biblio/532973 for MA956 ODS Steel Thermal Conductivity Data

clc;
clear all;
close all;


% ########################### Set Script Inputs
% UN Fuel Pellet Radius
r1 = 7.9; % mm

% MA 956 ODS Steel Radius
r2 = 8.4; % mm

% Maximum Fuel Centerline Temperature
Tm = 1600; % K

% Sodium Saturated Vapour Temperature Range
T_b = 800:100:1500; % K

% Initialise Fuel Linear Power Rating
q = NaN(length(T_b), 1); % W / m


% ########################### Import and Define Material Thermal Data
% MA956 ODS Steel Thermal Conductivity
k_MA956 = @ (T) ((1 ./ 60) .* (T - 273) + 10);


% ########################### Define Lambda Functions
% Define Helper Function to Solve for Fuel Pellet Surface Temperature
% This Equation is Derived by Integrating the Equation for UN Thermal Conductivity
T_s_UN_func = @ (x) ((Tm .^ 1.39...
                      - (1.39 ./ 1.41) .* (x ./ (4 .* pi))) .^ (1 ./ 1.39));

% Define Helper Function to Solve for Fuel Pellet Surface Temperature
% This Equation is Derived from Radial Heat Conduction Equation
T_s_TH_func = @ (x, T_b, G) (T_b + x .* (log(r2 / r1) ./ (2 .* pi .* k_MA956(T_b))));


% ########################### Calculate Linear Power Values
Tol = 1E-12;
for i = 1:length(T_b)
  t_b = T_b(i);

  try
    % Find Fuel Pin Linear Power Rating
    % Raise Function's Output to 1.39 to Avoid Imaginary Values
    tmp_func = @ (x) (real(T_s_UN_func(x) .^ 1.39) - (T_s_TH_func(x, t_b) .^ 1.39));
    [q_dash, val, info, out] = fzero(tmp_func, [100 1E6], optimset('TolX', Tol));
  catch
    continue;
  end_try_catch

  % Handle Case for Singular Point Convergence by Rejecting the q_dash Value
  if (info == -5)
    continue;
  endif

  q(i) = q_dash;
  clc;
  printf('Linear Power Calculation Progress >> %2.1f %%\n', ((100 * i) / length(T_b)));
endfor
clear Tol i j t_b q_dash tmp_func val info out;
clc;


% ########################### Plot Solution
% Plot for Reactor Linear Power with Fixed Range for Comparison
figure(1, 'name', 'Plot of Reactor Linear Power');
plot(T_b, q .* 1E-2, 'linewidth', 2);
hold on;
grid on;

##title(sprintf('Reactor Maximum Linear Power vs Sodium Vapour Temperature for Fuel Centerline Temperature of %.0f K', Tm));
xlabel('Saturated Sodium Vapour Temperature (K)');
ylabel('Reactor Maximum Linear Power (W/cm)');
set(gca, 'fontsize', 12);

% ########################### Clear Unnecessary Variables
clear -x r1 r2 Tm T_b q;

