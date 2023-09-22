% Script to Generate Heat Transfer Coefficient Tables for NSP Systems
% Refer ntrs.nasa.gov./citations./20060056311 for HeXe Coolant Data
% Refer doi.org/10.1016/0009-2614(75)80286-7 for He Coolant Data
% Refer physics.nist.gov./cgi-bin./Compositions./stand_alone.pl for Molar Masses

clc;
clear all;
close all;


% ########################### Set Script Inputs
G = 100:100:1500; % Kg / (m^2 * s)
T = 400:100:2000; % K
HeXe_MM = 40; % g / mol
D_h = 0.26354E-2; % m

% Fix Matrix Sizes for Calculation Vectorisation
G = ones(size(T, 2), 1) * G;
T  = T' * ones(1, size(G, 2));
HeXe_MM = HeXe_MM .* ones(size(G, 1), size(G, 2));
D_h = D_h .* ones(size(T, 1), size(T, 2));


% ########################### Import Processed HeXe Data from Johnson-2006
% Thermal Conductivity vs Temperature at Input Molar Mass
HeXe_k_Data = csvread('HeXe_k.csv');
HeXe_k = @ (T) (interp2(HeXe_k_Data(1, 2:end),...
               HeXe_k_Data(2:end, 1),...
               HeXe_k_Data(2:end, 2:end),...
               HeXe_MM, T));

% Viscosity vs Temperature at Input Molar Mass
HeXe_mu_Data = csvread('HeXe_mu.csv');
HeXe_mu = @ (T) (interp2(HeXe_mu_Data(1, 2:end),...
                HeXe_mu_Data(2:end, 1),...
                HeXe_mu_Data(2:end, 2:end),...
                HeXe_MM, T));

% Prandtl Number vs Temperature at Input Molar Mass
HeXe_Pr_Data = csvread('HeXe_Pr.csv');
HeXe_Pr = @ (T) (interp2(HeXe_Pr_Data(1, 2:end),...
                HeXe_Pr_Data(2:end, 1),...
                HeXe_Pr_Data(2:end, 2:end),...
                HeXe_MM, T));


% ########################### Import He Data from Jain-1975
% Thermal Conductivity vs Temperature
He_k_Data = csvread('He_k.csv');
He_k = @ (T) (interp1(He_k_Data(:, 1),...
             He_k_Data(:, 2),...
             T, 'spline'));

% Viscosity vs Temperature
He_mu_Data = csvread('He_mu.csv');
He_mu = @ (T) (interp1(He_mu_Data(:, 1),...
              He_mu_Data(:, 2),...
              T, 'spline'));

% Prandtl Number vs Temperature
He_Pr_Data = csvread('He_Pr.csv');
He_Pr = @ (T) (interp1(He_Pr_Data(:, 1),...
              He_Pr_Data(:, 2),...
              T, 'spline'));


% ########################### Calculate HTC Using Taylor-1988 for HeXe Mixtures
% Calculate Flow Reynolds Number
HeXe_Re = (G .* D_h ./ HeXe_mu(T));

% Enforce Reynolds Number Limits of HTC Correlation in Taylor-1988
HeXe_Re((HeXe_Re < 1.8E4) | (HeXe_Re > 6E4)) = NaN;

% Calculate HTC Values with T_W/T_b equal to 1
HeXe_HTC = (HeXe_k(T) ./ D_h) .* 0.023 .* (HeXe_Re .^ 0.8) .* (HeXe_Pr(T) .^ 0.65);


% ########################### Calculate HTC Using Taylor-1988 for He
% Calculate Flow Reynolds Number
He_Re = (G .* D_h ./ He_mu(T));

% Enforce Reynolds Number Limits of HTC Correlation in Taylor-1988
He_Re((He_Re < 1.8E4) | (He_Re > 6E4)) = NaN;

% Calculate HTC Values with T_W/T_b equal to 1
He_HTC = (He_k(T) ./ D_h) .* 0.023 .* (He_Re .^ 0.8) .* (He_Pr(T) .^ 0.65);


% ########################### Calculate Relative Difference Between HTCs
HTC_Diff = ((HeXe_HTC - He_HTC) ./ He_HTC) .* 100;


% ########################### Assemble and Save CSV Files
% Create CSV Data Container
CSV_Container = zeros(size(T, 1) + 1, size(G, 2) + 1);
CSV_Container(:, 1) = [-D_h(1, 1); T(:, 1)];
CSV_Container(1, :) = [-D_h(1, 1), G(1, :)];

% Save HeXe HTC Data to CSV
CSV_Container(2:end, 2:end) = HeXe_HTC; % W / m^2
csvwrite('HeXe_HTC.csv', CSV_Container);

% Save He HTC Data to CSV
CSV_Container(2:end, 2:end) = He_HTC; % W / m^2
csvwrite('He_HTC.csv', CSV_Container);


% ########################### Clear Unnecessary Variables
clear -x G T HeXe_MM D_h HeXe_HTC He_HTC HTC_Diff;

