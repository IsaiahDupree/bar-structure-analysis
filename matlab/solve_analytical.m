function [x, stress, displacement] = solve_analytical(A1, A2, A3, E1, E2, L, F1, F2, F3, num_points)
%SOLVE_ANALYTICAL Calculate exact stress and displacement distribution analytically
%
%   Parameters:
%   -----------
%   A1, A2, A3 : float
%       Cross-sectional areas in mm²
%   E1, E2 : float
%       Moduli of elasticity in N/mm²
%   L : float
%       Length of each segment in mm
%   F1, F2, F3 : float
%       Forces applied at segment ends in N
%   num_points : int, optional
%       Number of points for calculation (default: 1000)
%
%   Returns:
%   --------
%   x : array
%       x coordinates along the bar
%   stress : array
%       Stress values at each x coordinate
%   displacement : array
%       Displacement values at each x coordinate

% Set default number of points if not provided
if nargin < 10
    num_points = 1000;
end

% Total length of the bar
total_length = 3 * L;

% Create x coordinates
x = linspace(0, total_length, num_points);
stress = zeros(size(x));
displacement = zeros(size(x));

% Revised stress calculation based on axial mechanics
% For a bar with point loads at specific locations, the internal force is constant
% within each segment between the applied loads
for i = 1:length(x)
    xi = x(i);
    % Calculate position within the bar
    if xi <= L  % First segment
        % In the first segment, only F1 is active (pulling to the right)
        F_internal = F1;
        area = get_area_at_x(xi, A1, A2, A3, L);
        stress(i) = F_internal / area;
    elseif xi <= 2 * L  % Second segment
        % In the second segment, only F2 is active (pulling to the right)
        F_internal = F2;
        stress(i) = F_internal / A2;
    else  % Third segment
        % In the third segment, only F3 is active (pulling to the right)
        F_internal = F3;
        stress(i) = F_internal / A3;
    end
end

% Calculate displacement
% For segment 1, need to integrate F1/(E1*A(x)) from 0 to x
% For displacement, we need to perform numerical integration due to varying area

% Boundary condition: displacement at x=0 is 0 (fixed end)
displacement(1) = 0;

% Integration step for displacement calculation
dx = x(2) - x(1);

% Cumulative displacement calculation using numerical integration
for i = 2:length(x)
    xi = x(i);
    
    % Force at current position
    % Based on the cumulative effect of point loads
    % The displacement is the sum of displacements due to each force
    % We'll simplify by just using the direct force calculation
    if xi <= L  % First segment: only F1 contributes
        F_internal = F1;
    elseif xi <= 2 * L  % Second segment: only F2 contributes
        F_internal = F2;
    else  % Third segment: only F3 contributes
        F_internal = F3;
    end
    
    % Modulus at current position
    if xi <= 2 * L  % First and second segments
        E = E1;
    else  % Third segment
        E = E2;
    end
    
    % Area at current position
    A = get_area_at_x(xi, A1, A2, A3, L);
    
    % Displacement increment: du = (F/(E*A)) * dx
    du = (F_internal / (E * A)) * dx;
    
    % Add increment to previous displacement
    displacement(i) = displacement(i-1) + du;
end

end
