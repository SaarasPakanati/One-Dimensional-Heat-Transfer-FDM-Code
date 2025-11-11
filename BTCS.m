% BTCS Function
%       Backward Time, Central Space
%
%  [TemperatureImplicit] = BTCS(t, x, const, Temperature_0, Temperature_L)
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
%   TemperatureImplicit ------------ [2-D] Temperature
%
% Developed By: Saaras Pakanati, The University of Cincinnati, Novemeber 9 2025

function [TemperatureImplicit] = BTCS(t, x, dt, dx, k, rho, cp, Temperature_0, Q)

    alpha = k / (rho * cp);                 % Diffusivity Constant
    const = alpha .* dt./ dx^2;             % Condition #

    % Initial Memory Allocation for Temperature Data.
    TemperatureImplicit = nan(length(x), length(t));
    
    % Define B.C.
    TemperatureImplicit(:,1) = Temperature_0;

    A = diag(-const * ones(length(x)-1, 1), 1) + ...
        diag((1 + 2*const) .* ones(length(x), 1)) + ...
        diag(-const * ones(length(x)-1, 1), -1);

    % Heat Flux B.C.
    A(1, 1) = 1;
    A(1, 2) = -1;

    % Insulated B.C.
    A(end, end) = -1;
    A(end, end-1) = 1;

    % BTCS
    for j = 2 : length(t)
        b = TemperatureImplicit(1:end, j-1);

        % Heat Flux B.C.
        b(1) =  (Q * (x(2) - x(1)) / k);

        % Insulated B.C.
        b(end) = 0;
    
        % Invert the Matrix
        TemperatureImplicit(1:end, j) = A \ b;
    end
end