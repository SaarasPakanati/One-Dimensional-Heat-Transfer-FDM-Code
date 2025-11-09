%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       MATH 5106 Final Project                           %
%                                                                         %
% This is the main directory for the final project code developed in      %
% partial fulfillment of the requirements for MECH 5106                   %
%                                                                         %
% Function(s) Used:                                                       %
%   FTCS ------------------------- Forward Time Central Space FDM Solver  %
%   BTCS ------------------------- Backward Time Central Space FDM Solver %
%                                                                         %
% Developed By: Saaras Pakanati, The University of Cincinnati, 11/3/25    %
% Under the guidance of: Dr. Hyunjoong Kim, PhD                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initializing
clc; clear; close all;

% Physical Parameters
L = 1;                                  % [m] Length of the rod
t_tot = 1E3;                            % [s] Total time of simulation
k = 1E-5;                               % [W/m-K] Thermal Conductivity
rho = 1;                                % [kg/m^3] Density
cp = 1;                                 % [J/kg-K] Heat Capacity
Temperature_0 = 273;                    % [K] B.C Temperature at Node 1
Temperature_L = 373;                    % [K] B.C Temperature at Last Node
DeltaT = 25;                            % [K] Uncertainity in Temperature
DeltaA = 1E-6;                          % [] Uncertainity in Diffusivity Constant

% Numerical Parameters
dx = 1E-1;                              % [m] Spatial Discrtization
dt = 1E+1;                              % [s] Temporal Discretization
nRuns = 1E3;                            % [#] No. of Monte Carlo runs

% Array Initializations
x = 0 : dx : L;                         % X domain  
t = 0 : dt : t_tot;                     % T domain  
alpha = k / (rho * cp);                 % Diffusivity Constant
const = alpha .* dt./ dx^2;             % Condition #

% FTCS Method
[Temperature] = FTCS(t, x, const, Temperature_0, Temperature_L);

% BTCS Method
[TemperatureImplicit] = BTCS(t, x, const, Temperature_0, Temperature_L);

% Monte Carlo Method
[TemperatureMC] = MonteCarlo(dt, dx, t, x, alpha, Temperature_0, Temperature_L, nRuns, DeltaT, DeltaA);


plot(x, Temperature(:, end), 'r-', Linewidth = 2, DisplayName='FTCS'); hold on;
plot(x, TemperatureImplicit(:, end), 'k--', Linewidth = 2, DisplayName='BTCS');
plot(x, mean(TemperatureMC, 1), 'b:', 'LineWidth', 2, DisplayName='Monte-Carlo');
plot(x, mean(TemperatureMC, 1) + std(TemperatureMC, 0, 1), 'g--', 'LineWidth', 2, DisplayName='Error');
plot(x, mean(TemperatureMC, 1) - std(TemperatureMC, 0, 1), 'g--', 'LineWidth', 2, HandleVisibility='off');
legend(Location="best");
title(sprintf("t = %.2f s", t(end)), Interpreter="latex");
set(gca, 'FontWeight','bold','fontsize',14, 'FontName', 'Serif');
xlabel('x [$m$]', Interpreter='latex');
ylabel('T [$K$]', Interpreter='latex');
grid on;