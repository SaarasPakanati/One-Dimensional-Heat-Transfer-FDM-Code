% MonteCarlo Function
%       Monte-Carlo Method to solve for Temperature
%
%  [TemperatureMC, maxTemp] = MonteCarlo(t, x, const, Temperature_0, Q, nRuns)
%
% Input(s):
%   t ------------------------------ [s] Temporal Coordinate
%   x ------------------------------ [m] Spatial Coordinate
%   dt ----------------------------- [s] Temporal Discretization
%   dx ----------------------------- [m] Spatial Discretization
%   k ------------------------------ [W/m-K] Thermal Conductivity
%   rho ---------------------------- [kg/m^3] Density
%   cp ----------------------------- [J/kg-K] Heat Capacity
%   Temperature_0 ------------------ [K] Value at Node 0, at t = 0
%   Q ------------------------------ [W/m^2] Heat Flux
%   nRuns -------------------------- [#] Number of MonteCarlo runs
%   T_crit ------------------------- [K] Critical Temperature
%
% Output(s):
%   dataTemp ----------------------- [2-D] Point of Interests in Temperature Data
%
% Developed By: Saaras Pakanati, The University of Cincinnati, Novemeber 3 2025

function [dataTemp] = MonteCarlo(t, x, dt, dx, k, rho, cp, Temperature_0, Q, nRuns, T_crit)

    maxTemp = [nan, nan];

    for l = 1:nRuns

        k_temp = normrnd(k, 0.05);
        
        while k_temp < 0
            k_temp = normrnd(k, 0.05);
        end

        [TemperatureMC] = BTCS(t, x, dt, dx, k_temp, rho, cp, Temperature_0, Q);

        maxTemp = [maxTemp; [k_temp, max(TemperatureMC(:, end))]];

        disp(l*100/nRuns);
    end

    dataTemp = [maxTemp(:, 1), maxTemp(:, 2)];


end