% MonteCarlo Function
%       Monte-Carlo Method to solve for Temperature
%
%  [TemperatureMC] = MonteCarlo(t, x, const, Temperature_0, Temperature_L)
%
% Input(s):
%   dt ----------------------------- [s] Discretized Temporal Coordinate
%   dx ----------------------------- [m] Discretized Spatial Coordinate
%   t ------------------------------ [s] Temporal Coordinate
%   x ------------------------------ [m] Spatial Coordinate
%   alpha -------------------------- [] Diffusive Constant
%   Temperature_0 ------------------ [K] Value at Node 0
%   Temperature_L ------------------ [K] Value at Node end
%   nRuns -------------------------- [#] No. of Monte Carlo Runs
%   DeltaT ------------------------- [K] Variation in Temperature
%   DeltaA ------------------------- [] Variation in Diffusivity Constant
%
% Output(s):
%   TemperatureMC -------------------- [2-D] Temperature
%
% Developed By: Saaras Pakanati, The University of Cincinnati, Novemeber 3 2025

function [TemperatureMC] = MonteCarlo(dt, dx, t, x, alpha, Temperature_0, Temperature_L, nRuns, DeltaT, DeltaA)

    % Initial Memory Allocation for Temperature Data.
    TemperatureMC = zeros(nRuns, length(x));

    for k = 1:nRuns
        % Initialize Physical Parameter Values
        a =  alpha         + DeltaA * randn();
        T0 = Temperature_0 + DeltaT * randn();
        TL = Temperature_L + DeltaT * randn();
        
        % Initialize Numerical Parameter Values
        const = a * dt / dx^2;
        
        % FTCS Solve
        [TemperatureTemp] = FTCS(t, x, const, T0, TL);

        % Save Data from Temporary File.
        TemperatureMC(k,:) = TemperatureTemp(:, end);
    end

end