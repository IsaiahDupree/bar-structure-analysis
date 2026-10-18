function displayAnalyticalSolution(P)
% DISPLAYANALYTICALSOLUTION Display the analytical solution for the bar structure
%
% Parameters:
% P - Structure containing problem parameters

fprintf('\n========== ANALYTICAL SOLUTION ==========\n\n');

% Calculate total reaction force at fixed end
R = sum(P.F);
fprintf('Total reaction at fixed end: R = %.1f kN\n', R/1000);

% Internal forces in each segment
N1 = R;                    % Segment 1: N1 = R
N2 = R - P.F(1);           % Segment 2: N2 = R - F1
N3 = R - P.F(1) - P.F(2);  % Segment 3: N3 = R - F1 - F2

fprintf('Internal forces:\n');
fprintf('  Segment 1: N1 = %.1f kN\n', N1/1000);
fprintf('  Segment 2: N2 = %.1f kN\n', N2/1000);
fprintf('  Segment 3: N3 = %.1f kN\n', N3/1000);

% Stress calculations
sigma1_start = N1 / P.A(1);
sigma1_end = N1 / P.A(2);
sigma2 = N2 / P.A(2);
sigma3 = N3 / P.A(3);

fprintf('\nStress distribution:\n');
fprintf('  Segment 1: σ₁ varies from %.1f MPa to %.1f MPa\n', sigma1_start, sigma1_end);
fprintf('  Segment 2: σ₂ = %.1f MPa (constant)\n', sigma2);
fprintf('  Segment 3: σ₃ = %.1f MPa (constant)\n', sigma3);

% Displacement calculations
% For segment 1 with linearly varying cross-section
% The formula is: δ₁ = (N₁*L)/(E₁(A₂-A₁))*ln(A₂/A₁)
delta1 = (N1 * P.L) / (P.E(1) * (P.A(2) - P.A(1))) * log(P.A(2) / P.A(1));

% For segments 2 and 3 with constant cross-section
% The formula is: δᵢ = (Nᵢ*L)/(Eᵢ*Aᵢ)
delta2 = (N2 * P.L) / (P.E(1) * P.A(2));
delta3 = (N3 * P.L) / (P.E(2) * P.A(3));

% Total end displacement
total_displacement = delta1 + delta2 + delta3;

fprintf('\nAxial displacements:\n');
fprintf('  Segment 1: δ₁ = %.2f mm\n', delta1);
fprintf('  Segment 2: δ₂ = %.2f mm\n', delta2);
fprintf('  Segment 3: δ₃ = %.2f mm\n', delta3);
fprintf('  Total end displacement: %.2f mm\n\n', total_displacement);

end
