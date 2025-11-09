% BTCS Function
%       Backward Time, Central Space
%
%  [TemperatureImplicit] = BTCS(t, x, const, Temperature_0, Temperature_L)
%
% Input(s):
%   t ------------------------------ [s] Temporal Coordinate
%   x ------------------------------ [m] Spatial Coordinate
%   const -------------------------- [] Diffusive Constant
%   Temperature_0 ------------------ [K] Value at Node 0
%   Temperature_L ------------------ [K] Value at Node end
%
% Output(s):
%   TemperatureImplicit ------------ [2-D] Temperature
%
% Developed By: Saaras Pakanati, The University of Cincinnati, Novemeber 3 2025

function [TemperatureImplicit] = BTCS(t, x, const, Temperature_0, Temperature_L)

    % Initial Memory Allocation for Temperature Data.
    TemperatureImplicit = nan(length(x), length(t));
    
    % Define B.Cs.
    TemperatureImplicit(1,:) = Temperature_0;
    TemperatureImplicit(:,1) = Temperature_0;
    TemperatureImplicit(end,:) = Temperature_L;
    
    % Diagonal Matrix Definition
    A = diag(-const * ones(length(x)-3, 1), 1) + ...
        diag((1 + 2*const) .* ones(length(x)-2, 1)) + ...
        diag(-const * ones(length(x)-3, 1), -1);
    
    % BTCS
    for j = 2 : length(t)
        b = TemperatureImplicit(2:end-1, j-1);
    
        b(1) = b(1) + const * TemperatureImplicit(1, j-1);
        b(end) = b(end) + const * TemperatureImplicit(end, j-1);
    
        TemperatureImplicit(2:end-1, j) = A \ b;
    end
end