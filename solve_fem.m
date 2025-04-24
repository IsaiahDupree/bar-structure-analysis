function [x_nodal, nodal_displacements, element_stresses, error] = solve_fem(A1, A2, A3, E1, E2, L, F1, F2, F3, num_elements_per_segment)
%SOLVE_FEM Solve using finite element method with 2-node elements
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
%   num_elements_per_segment : int
%       Number of elements per segment
%
%   Returns:
%   --------
%   x_nodal : array
%       x coordinates of nodes
%   nodal_displacements : array
%       Displacement at each node
%   element_stresses : array
%       Stress in each element
%   error : array
%       Error compared to analytical solution (at element centers)

% Total length of the bar
total_length = 3 * L;

% Total number of elements
num_elements = 3 * num_elements_per_segment;

% Number of nodes (for 2-node elements)
num_nodes = num_elements + 1;

% Element length
element_length = total_length / num_elements;

% Create node coordinates
x_nodal = linspace(0, total_length, num_nodes);

% Initialize global stiffness matrix
K_global = zeros(num_nodes, num_nodes);

% Initialize force vector
F_global = zeros(num_nodes, 1);

% Apply point loads at the nodes corresponding to segment ends
% This matches our revised analytical model where forces are applied at segment junctions

% F1 at the end of segment 1 (node at x = L)
F_global(num_elements_per_segment + 1) = F1;

% F2 at the end of segment 2 (node at x = 2L)
F_global(2 * num_elements_per_segment + 1) = F2;

% F3 at the end of segment 3 (node at x = 3L)
F_global(end) = F3;

% List to store element stiffness matrices and properties
element_properties = cell(num_elements, 1);

% Assemble global stiffness matrix
for i = 1:num_elements
    % Element nodes
    node1 = i;
    node2 = i + 1;
    
    % Element length
    Le = element_length;
    
    % Calculate properties at element center
    x_center = (x_nodal(node1) + x_nodal(node2)) / 2;
    
    % Get modulus and area at element center
    if x_center <= 2 * L  % First and second segments
        E = E1;
    else  % Third segment
        E = E2;
    end
    
    if x_center <= L  % First segment
        % Linear interpolation from A1 to A2
        A = A1 + (A2 - A1) * (x_center / L);
    elseif x_center <= 2 * L  % Second segment
        A = A2;
    else  % Third segment
        A = A3;
    end
    
    % Element stiffness matrix: k = (A*E/L) * [1 -1; -1 1]
    k = (A * E / Le) * [1, -1; -1, 1];
    
    % Store element stiffness matrix and properties
    element_properties{i} = struct('k', k, 'node1', node1, 'node2', node2, 'A', A, 'E', E, 'Le', Le);
    
    % Assemble into global stiffness matrix
    K_global(node1, node1) = K_global(node1, node1) + k(1, 1);
    K_global(node1, node2) = K_global(node1, node2) + k(1, 2);
    K_global(node2, node1) = K_global(node2, node1) + k(2, 1);
    K_global(node2, node2) = K_global(node2, node2) + k(2, 2);
end

% Apply boundary condition (fixed at x=0)
% Remove first row and column from stiffness matrix
K_reduced = K_global(2:end, 2:end);
F_reduced = F_global(2:end);

% Solve for nodal displacements
% KU = F => U = K^-1 * F
U_reduced = K_reduced \ F_reduced;

% Add zero displacement at fixed end
nodal_displacements = zeros(num_nodes, 1);
nodal_displacements(2:end) = U_reduced;

% Calculate element stresses
element_stresses = zeros(num_elements, 1);

for i = 1:num_elements
    % Get element properties
    props = element_properties{i};
    node1 = props.node1;
    node2 = props.node2;
    A = props.A;
    E = props.E;
    Le = props.Le;
    
    % Element displacement
    u1 = nodal_displacements(node1);
    u2 = nodal_displacements(node2);
    
    % Strain = du/dx
    strain = (u2 - u1) / Le;
    
    % Stress = E * strain
    element_stresses(i) = E * strain;
end

% Calculate error compared to analytical solution
[x_analytical, stress_analytical, ~] = solve_analytical(A1, A2, A3, E1, E2, L, F1, F2, F3);

% Element center coordinates
x_element_centers = (x_nodal(1:end-1) + x_nodal(2:end)) / 2;

% Interpolate analytical solution at element centers
stress_analytical_at_centers = interp1(x_analytical, stress_analytical, x_element_centers);

% Calculate relative error with improved handling
error = zeros(size(element_stresses));
for i = 1:length(error)
    if abs(stress_analytical_at_centers(i)) > 1.0  % Only use relative error for significant stress values
        error(i) = (element_stresses(i) - stress_analytical_at_centers(i)) / (abs(stress_analytical_at_centers(i)) + 1e-10);
    else
        % Use absolute error if stress is very small
        error(i) = element_stresses(i) - stress_analytical_at_centers(i);
    end
end

end
