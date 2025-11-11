% FTCS Function
%       Forward Time, Central Space
%
%  [Temperature] = FTCS(t, x, const, Temperature_0, Temperature_L)
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
%
% Output(s):
%   Temperature -------------------- [2-D] Temperature
%
% Developed By: Saaras Pakanati, The University of Cincinnati, Novemeber 9 2025

function [Temperature] = FTCS(t, x, dt, dx, k, rho, cp, Temperature_0, Q)

    alpha = k / (rho * cp);                 % Diffusivity Constant
    const = alpha .* dt./ dx^2;             % Condition #

    % Initial Memory Allocation for Temperature Data.
    Temperature = nan(length(x), length(t));
    
    % Define B.C.
    Temperature(:, 1) = Temperature_0;
    
    % FTCS
    for j = 1 : length(t) - 1
        for i = 2 : length(x) - 1
    
            Temperature(i, j+1) = (const) .* (Temperature(i+1, j) ...
                - 2.* Temperature(i, j) + Temperature(i-1, j)) + Temperature(i, j);
    
        end

        % Heat Flux B.C.
        Temperature(1, j+1) = Temperature(2, j+1) + (Q * (x(2) - x(1)) / k);

        % Insulated B.C.
        Temperature(end, j+1) = Temperature(end-1, j+1);

        % plot(x, Temperature(:, j+1), 'r-', Linewidth = 2, DisplayName='FTCS');
        % title(sprintf("t = %.2f s", t(j)), Interpreter="latex");drawnow();

    end
end