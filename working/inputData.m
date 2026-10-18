function P = inputData
% INPUTDATA Define parameters for the bar structure problem
%
% Returns:
% P - Structure containing problem parameters

% Geometry and material parameters
P.L   = 500;                    % mm
P.F   = [20e3 40e3 20e3];       % N
P.A   = [200 100 50];           % mm^2
P.E   = [130e3 200e3];          % MPa
P.nEl = [4 4 4];                % starting mesh (can refine)

end
