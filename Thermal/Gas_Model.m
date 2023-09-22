% Script to Find Bounding Fuel Element Linear Power for Gas Cooled NSP Systems
% Refer ntrs.nasa.gov./citations./20060056311 for HeXe Coolant Data
% Refer doi.org/10.1016/0009-2614(75)80286-7 for He Coolant Data
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

% Bulk Coolant Temperature Range
T_b = 800:100:1400; % K

% Coolant Mass Flux Range
G = 100:100:1000; % Kg / (m^2 * s)

% Initialise Fuel Linear Power Rating
q = NaN(length(T_b), length(G)); % W / m


% ########################### Import and Define Material Thermal Data
% Coolant HTC Values from Johnson-2006 and Taylor-1988
HTC_Data = csvread('HeXe_HTC.csv');
% Replace NaN Values with a Dummy Number for interp2
HTC_Data(isnan(HTC_Data)) = 1;
HTC = @ (G, T) (interp2(HTC_Data(1, 2:end),...
                        HTC_Data(2:end, 1),...
                        HTC_Data(2:end, 2:end),...
                        G, T));

##% Coolant HTC Values from Jain-1975 and Taylor-1988
##HTC_Data = csvread('He_HTC.csv');
##% Replace NaN Values with a Dummy Number for interp2
##HTC_Data(isnan(HTC_Data)) = 1;
##HTC = @ (G, T) (interp2(HTC_Data(1, 2:end),...
##                        HTC_Data(2:end, 1),...
##                        HTC_Data(2:end, 2:end),...
##                        G, T));

% MA956 ODS Steel Thermal Conductivity
k_MA956 = @ (T) ((1 ./ 60) .* (T - 273) + 10);


% ########################### Define Lambda Functions
% Define Helper Function to Solve for Fuel Pellet Surface Temperature
% This Equation is Derived by Integrating the Equation for UN Thermal Conductivity
T_s_UN_func = @ (x) ((Tm .^ 1.39...
                      - (1.39 ./ 1.41) .* (x ./ (4 .* pi))) .^ (1 ./ 1.39));

% Define Helper Function to Solve for Fuel Pellet Surface Temperature
% This Equation is Derived from Radial Heat Conduction and Convection Equations
T_s_TH_func = @ (x, T_b, G) (T_b + x .* (...
                             (log(r2 / r1) ./ (2 .* pi .* k_MA956(T_b)))...
                             + (1 ./ (2 .* pi .* r2 .* 1E-3 .* HTC(G, T_b)))));


% ########################### Calculate Linear Power Values
Tol = 1E-12;
for i = 1:length(T_b)
  t_b = T_b(i);
  for j = 1:length(G)
    g = G(j);

    try
      % Find Fuel Pin Linear Power Rating
      % Raise Function's Output to 1.39 to Avoid Imaginary Values
      tmp_func = @ (x) (real(T_s_UN_func(x) .^ 1.39) - (T_s_TH_func(x, t_b, g) .^ 1.39));
      [q_dash, val, info, out] = fzero(tmp_func, [100 1E6], optimset('TolX', Tol));
    catch
      continue;
    end_try_catch

    % Handle Case for Singular Point Convergence by Rejecting the q_dash Value
    if (info == -5)
      continue;
    endif

    q(i, j) = q_dash;
  endfor
  clc;
  printf('Linear Power Calculation Progress >> %2.1f %%\n', ((100 * i) / length(T_b)));
endfor
clear Tol i j t_b g q_dash tmp_func val info out;
clc;


% ########################### Plot Solution
% Plot for Reactor Linear Power with Fixed Range for Comparison
figure(1, 'name', 'Reactor Linear Power');
% Define Linestyles and Marker Types to Cycle to Through
mrkrs = ['+', 'o', '*', '.', 'x', 's', 'd', '^', 'v', '>', '<'];
mrkrs_cycler = @ (i) (mrkrs(mod(i, 11) + 1));
lines = ['-', '--', ':', '-.'];
line_cycler = @ (i) (lines((mod(i, 4) + 1)));
for i = 1:length(T_b)
  plot(G, q(i, :) .* 1E-2, sprintf('%s%s;HeXe Coolant Temperature: %.0f K;',...
       mrkrs_cycler(i), line_cycler(i), T_b(i)),...
       'linewidth', 2);
  hold on;
  grid on;
endfor

##title(sprintf('Reactor Maximum Linear Power vs Core Mass Flux for Fuel Centerline Temperature of %.0f K', Tm));
xlabel('Core Mass Flux (Kg/m^2s)');
ylabel('Reactor Maximum Linear Power (W/cm)');
ylim([0, 700]);
legend('location', 'southwest');
set(gca, 'fontsize', 12);

% ########################### Clear Unnecessary Variables
clear -x r1 r2 Tm T_b G q HTC_Data;

