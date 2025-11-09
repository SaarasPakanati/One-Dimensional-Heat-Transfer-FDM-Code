% FTCS Function
%       Forward Time, Central Space
%
%  [Temperature] = FTCS(t, x, const, Temperature_0, Temperature_L)
%
% Input(s):
%   t ------------------------------ [s] Temporal Coordinate
%   x ------------------------------ [m] Spatial Coordinate
%   const -------------------------- [] Diffusive Constant
%   Temperature_0 ------------------ [K] Value at Node 0
%   Temperature_L ------------------ [K] Value at Node end
%
% Output(s):
%   Temperature -------------------- [2-D] Temperature
%
% Developed By: Saaras Pakanati, The University of Cincinnati, Novemeber 3 2025

function [Temperature] = FTCS(t, x, const, Temperature_0, Temperature_L)

    % Initial Memory Allocation for Temperature Data.
    Temperature = nan(length(x), length(t));
    
    % Define B.Cs.
    Temperature(1,:) = Temperature_0;
    Temperature(:,1) = Temperature_0;
    Temperature(end,:) = Temperature_L;
    
    % FTCS
    for j = 1 : length(t) - 1
        for i = 2 : length(x) - 1
    
            Temperature(i, j+1) = (const) .* (Temperature(i+1, j) ...
                - 2.* Temperature(i, j) + Temperature(i-1, j)) + Temperature(i, j);
    
        end
    end
end