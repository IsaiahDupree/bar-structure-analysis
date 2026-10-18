function [mesh, areas, Es] = generateMesh(P)
% GENERATEMESH Generate mesh for the composite bar
%
% Parameters:
% P - Structure containing problem parameters
%
% Returns:
% mesh - Matrix with columns [node1, node2, x0, length]
% areas - Vector of element cross-sectional areas
% Es - Vector of element elastic moduli

% Initialize arrays
node = 0; 
mesh = []; 
areas = []; 
Es = [];

% Extract cross-sectional areas
A1 = P.A(1); 
A2 = P.A(2); 
A3 = P.A(3);

% Generate mesh for each segment
for s = 1:3
    n = P.nEl(s);  % Number of elements in this segment
    x0 = (s-1)*P.L;  % Starting x-coordinate for this segment
    le = P.L/n;  % Element length
    
    for e = 1:n
        % Add element: [node1, node2, x0, length]
        mesh(end+1,:) = [node+e-1, node+e, x0+(e-1)*le, le];
        
        if s == 1  % First segment - linear taper
            xmid = x0 + (e-0.5)*le;  % Midpoint of element
            areas(end+1) = A1 + (A2-A1)*(xmid/P.L);
            Es(end+1) = P.E(1);
        elseif s == 2  % Second segment - constant A2
            areas(end+1) = A2;
            Es(end+1) = P.E(1);
        else  % Third segment - constant A3
            areas(end+1) = A3;
            Es(end+1) = P.E(2);
        end
    end
    
    % Update node counter for next segment
    node = node + n;
end

end
