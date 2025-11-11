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
L = 5E-2;                               % [m] Length of the cross-section
t_tot = 3E+2;                           % [s] Total time of simulation
k = 5E-1;                               % [W/m-K] Thermal Conductivity
Q = 1E+5;                               % [W/m^2] Heat Flux
rho = 3E+2;                             % [kg/m^3] Density
cp = 1E+3;                              % [J/kg-K] Heat Capacity
Temperature_0 = 300;                    % [K] B.C Temperature, Initially
T_crit = 7E+3;                          % [K] Critical Temperature, above which material will fail
DeltaT = 25;                            % [K] Uncertainity in Temperature
DeltaA = 1E-6;                          % [] Uncertainity in Diffusivity Constant

% Numerical Parameters
dx = 5E-3;                              % [m] Spatial Discrtization
dt = 2E-0;                              % [s] Temporal Discretization
nRuns = 1E3;                            % [#] No. of Monte Carlo runs

% Array Initializations
x = 0 : dx : L;                         % X domain  
t = 0 : dt : t_tot;                     % T domain  
alpha = k / (rho * cp);                 % Diffusivity Constant
const = alpha .* dt./ dx^2;             % Condition #

% Check Convergence

fprintf("The Diffusion coefficient = %.2f \n", const);

exitPrompt = 'N';
exitPromptMC = 'N';

if const > 1
    fprintf("Unstable solution expected because %.2f > 1 \n", const);
    exitPrompt = input("Would you like to exit the program? [Y/N]:  ");
else
    fprintf("Stable solution expected because %.2f < 1 \n", const);
end

if exitPrompt == 'N'
    % FTCS Method
    [Temperature] = FTCS(t, x, dt, dx, k, rho, cp, Temperature_0, Q);
    
    % BTCS Method
    [TemperatureImplicit] = BTCS(t, x, dt, dx, k, rho, cp, Temperature_0, Q);

    % Plot
    figure();
    plot(x, Temperature(:, end), 'r-', Linewidth = 2, DisplayName='FTCS'); hold on;
    plot(x, TemperatureImplicit(:, end), 'k--', Linewidth = 2, DisplayName='BTCS');
    legend(Location="best");
    title(sprintf("t = %.2f s | %.5f", t(end), const), Interpreter="latex");
    set(gca, 'FontWeight','bold','fontsize',14, 'FontName', 'Serif');
    xlabel('x [$m$]', Interpreter='latex');
    ylabel('T [$K$]', Interpreter='latex');
    grid on;
    
    exitPromptMC = input("Would you like to proceed to the Monte Carlo Simulation? [Y/N]:  ");
    % exitPromptMC = 'Y';
end

if exitPromptMC == 'Y'
    % Monte Carlo Method
    [dataTemp] = MonteCarlo(t, x, dt, dx, k, rho, cp, Temperature_0, Q, nRuns, T_crit);

    % Plot Histogram
    figure;
    histogram(dataTemp(:, 2), 200, 'Normalization', 'probability', 'FaceColor', 'b', 'FaceAlpha', 0.7, HandleVisibility='off');
    hold on;
    xline(T_crit, 'r--', 'LineWidth', 2, 'Label', 'Critical Temperature', HandleVisibility='off');
    xline(mean(dataTemp(:, 2), 'omitnan'), 'k-', 'LineWidth', 2, 'Label', 'Mean', HandleVisibility='off');
    xlabel('Maximum Temperature [K]');
    ylabel('Probability [-]');
    title('Temperature Distribution from Monte Carlo Simulation');
    set(gca, 'FontWeight','bold','fontsize',14, 'FontName', 'Serif');
    grid on;
    
    % Plot CDF
    figure;
    [f, x] = ecdf(dataTemp(:, 2));
    plot(x, f*100, 'b-', 'LineWidth', 2, HandleVisibility='off');
    hold on;
    xline(T_crit, 'r--', 'LineWidth', 2, 'Label', 'Critical Temperature', 'LabelVerticalAlignment', 'bottom', HandleVisibility='off');
    yline(68.3, 'g--', 'LineWidth', 1.5, 'Label', '1 \sigma', DisplayName='1 \sigma');
    yline(95.45, 'g--', 'LineWidth', 1.5, 'Label', '2 \sigma', DisplayName='2 \sigma');
    yline(99.7, 'g--', 'LineWidth', 1.5, 'Label', '3 \sigma', DisplayName='3 \sigma');
    xlabel('Temperature [K]');
    ylabel('Cumulative Probability [%]');
    title('Cummulative Density Function');
    set(gca, 'FontWeight','bold','fontsize',14, 'FontName', 'Serif');
    grid on;
    axis padded;
    legend(Location="south");
    
    % Box Plot
    figure;
    boxplot(dataTemp(:, 2), 'Labels', {'Max Temperature'});
    hold on;
    yline(200, 'r--', 'LineWidth', 2);
    ylabel('Temperature [K]');
    title('Box Plot');
    set(gca, 'FontWeight','bold','fontsize',14, 'FontName', 'Serif');
    grid on;
    
    % Thermal Conductivity v/s Max Temperature Plot
    figure;
    scatter(dataTemp(:, 1), dataTemp(:, 2), 20, 'filled', 'MarkerFaceAlpha', 0.5);hold on;
    yline(T_crit, 'r--', 'LineWidth', 2, 'Label', 'Critical Temperature');
    xlabel('Thermal Conductivity k [W/m·K]');
    ylabel('Max Temperature [K]');
    title('Sensitivity: k v/s T_{max}');
    grid on;
    set(gca, 'FontWeight','bold','fontsize',14, 'FontName', 'Serif');

    % Plot Probability of Temperature over Critical Temperature
    figure();
    thresholds = linspace(min(dataTemp(:, 2)), max(dataTemp(:, 2)), 100);
    exceedance_prob = zeros(size(thresholds));
    
    for i = 1:length(thresholds)
        exceedance_prob(i) = sum(dataTemp(:, 2) > thresholds(i)) / length(dataTemp(:, 2)) * 100;
    end
    
    plot(thresholds, exceedance_prob, 'b-', 'LineWidth', 2);
    hold on;
    xline(T_crit, 'r--', 'LineWidth', 2, 'Label', 'Critical Temperature');
    xlabel('Temperature Threshold [K]');
    ylabel('Exceedance Probability [%]');
    title('Risk Curve');
    grid on;
    set(gca, 'YScale', 'log', 'FontWeight','bold','fontsize',14, 'FontName', 'Serif');
    
end

%%
plot(x, Temperature(:, end), 'r-', Linewidth = 2, DisplayName='FTCS'); hold on;
plot(x, TemperatureImplicit(:, end), 'k--', Linewidth = 2, DisplayName='BTCS');
% plot(x, mean(TemperatureMC, 1), 'b:', 'LineWidth', 2, DisplayName='Monte-Carlo');
% plot(x, mean(TemperatureMC, 1) + std(TemperatureMC, 0, 1), 'g--', 'LineWidth', 2, DisplayName='Error');
% plot(x, mean(TemperatureMC, 1) - std(TemperatureMC, 0, 1), 'g--', 'LineWidth', 2, HandleVisibility='off');
legend(Location="best");
title(sprintf("t = %.2f s", t(end)), Interpreter="latex");
set(gca, 'FontWeight','bold','fontsize',14, 'FontName', 'Serif');
xlabel('x [$m$]', Interpreter='latex');
ylabel('T [$K$]', Interpreter='latex');
grid on;