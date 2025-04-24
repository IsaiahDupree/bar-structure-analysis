function [U, stress] = solveBar(mesh, A, E, P)
% SOLVEBAR Solve the bar structure using finite element method
%
% Parameters:
% mesh - Matrix with columns [node1, node2, x0, length]
% A - Vector of element cross-sectional areas
% E - Vector of element elastic moduli
% P - Structure containing problem parameters
%
% Returns:
% U - Vector of nodal displacements
% stress - Vector of element stresses

% Get number of nodes
nNode = max(mesh(:,2)) + 1;

% Initialize global stiffness matrix
K = zeros(nNode);

% Assemble global stiffness matrix
for e = 1:size(mesh, 1)
    % Element nodes (adding 1 because MATLAB is 1-indexed)
    n1 = mesh(e,1) + 1;
    n2 = mesh(e,2) + 1;
    
    % Element length
    le = mesh(e,4);
    
    % Element stiffness matrix: k = (A*E/L) * [1 -1; -1 1]
    ke = A(e) * E(e) / le * [1, -1; -1, 1];
    
    % Assemble into global stiffness matrix
    K([n1 n2], [n1 n2]) = K([n1 n2], [n1 n2]) + ke;
end

% Initialize global force vector
F = zeros(nNode, 1);

% Apply forces at segment junctions
F(end) = P.F(3);     % F3 at end of segment 3
F(end-1) = P.F(2);   % F2 at junction of segments 2 and 3
F(end-2) = P.F(1);   % F1 at junction of segments 1 and 2

% Apply boundary condition (fixed at x=0)
K(1,:) = 0; 
K(:,1) = 0; 
K(1,1) = 1;  
F(1) = 0;

% Solve for nodal displacements: KU = F
U = K \ F;

% Calculate element stresses
stress = zeros(size(mesh, 1), 1);

for e = 1:size(mesh, 1)
    % Element nodes
    n1 = mesh(e,1) + 1;
    n2 = mesh(e,2) + 1;
    
    % Element length
    le = mesh(e,4);
    
    % Stress = E * strain = E * (u2-u1)/L
    stress(e) = E(e) / le * (U(n2) - U(n1));
end

end
