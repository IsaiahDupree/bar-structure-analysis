%% Enhanced Bar Structure Analysis - Main Script
% This script analyzes a composite bar structure using both analytical and
% finite element methods, with enhanced features matching the Python implementation

clear all;
close all;
clc;

%% Define parameters
% Problem parameters
A1 = 200;  % mm²
A2 = 100;  % mm²
A3 = 50;   % mm²
E1 = 130;  % GPa
E2 = 200;  % GPa
L = 500;   % mm
F1 = 20;   % kN
F2 = 40;   % kN
F3 = 20;   % kN

% Convert units to consistent system (N and mm)
E1 = E1 * 1000;  % Convert from GPa to N/mm²
E2 = E2 * 1000;  % Convert from GPa to N/mm²
F1 = F1 * 1000;  % Convert from kN to N
F2 = F2 * 1000;  % Convert from kN to N
F3 = F3 * 1000;  % Convert from kN to N

%% Display header
fprintf('========== BAR STRUCTURE ANALYSIS ==========\n');
fprintf('Composite Bar with Variable Cross-Section\n');
fprintf('===========================================\n');
fprintf('Parameters:\n');
fprintf('A1 = %.0f mm², A2 = %.0f mm², A3 = %.0f mm²\n', A1, A2, A3);
fprintf('E1 = %.1f GPa, E2 = %.1f GPa\n', E1/1000, E2/1000);
fprintf('L = %.0f mm\n', L);
fprintf('F1 = %.1f kN, F2 = %.1f kN, F3 = %.1f kN\n\n', F1/1000, F2/1000, F3/1000);

%% 1. ANALYTICAL SOLUTION
fprintf('1. ANALYTICAL (CLOSED-FORM) SOLUTION:\n');
fprintf('%s\n', repmat('-', 1, 40));

% Calculate total reaction at fixed end
R = F1 + F2 + F3;
fprintf('Total reaction at fixed end: R = %.1f kN\n', R/1000);

% Internal forces in each segment
N1 = R;
N2 = R - F1;
N3 = R - F1 - F2;

fprintf('Internal forces:\n');
fprintf('  Segment 1: N1 = %.1f kN\n', N1/1000);
fprintf('  Segment 2: N2 = %.1f kN\n', N2/1000);
fprintf('  Segment 3: N3 = %.1f kN\n', N3/1000);

% Stress calculations
sigma1_start = N1 / A1;
sigma1_end = N1 / A2;
sigma2 = N2 / A2;
sigma3 = N3 / A3;

fprintf('Stress distribution:\n');
fprintf('  Segment 1: σ₁ varies from %.1f MPa to %.1f MPa\n', sigma1_start, sigma1_end);
fprintf('  Segment 2: σ₂ = %.1f MPa (constant)\n', sigma2);
fprintf('  Segment 3: σ₃ = %.1f MPa (constant)\n', sigma3);

% Displacement calculations
% For segment 1 with linearly varying cross-section
% The formula is: δ₁ = (N₁*L)/(E₁(A₂-A₁))*ln(A₂/A₁)
delta1 = (N1 * L) / (E1 * (A2 - A1)) * log(A2 / A1);

% For segments 2 and 3 with constant cross-section
% The formula is: δᵢ = (Nᵢ*L)/(Eᵢ*Aᵢ)
delta2 = (N2 * L) / (E1 * A2);
delta3 = (N3 * L) / (E2 * A3);

% Total end displacement
total_displacement = delta1 + delta2 + delta3;

fprintf('Axial displacements:\n');
fprintf('  Segment 1: δ₁ = %.2f mm\n', delta1);
fprintf('  Segment 2: δ₂ = %.2f mm\n', delta2);
fprintf('  Segment 3: δ₃ = %.2f mm\n', delta3);
fprintf('  Total end displacement: %.2f mm\n\n', total_displacement);

%% 2. FINITE ELEMENT SOLUTION
fprintf('2. FINITE ELEMENT ANALYSIS:\n');
fprintf('%s\n', repmat('-', 1, 40));

% Initial number of elements per segment
num_elements_per_segment = 4;

% Calculate analytical solution
tic;
fprintf('Computing analytical solution...\n');
[x_analytical, stress_analytical, displacement_analytical] = solve_analytical(A1, A2, A3, E1, E2, L, F1, F2, F3);
analytical_time = toc;
fprintf('Analytical solution computed in %.2f seconds\n', analytical_time);

%% Use finite element method with different element counts until error is below 5%
max_error = inf;
max_iterations = 5;
iterations = 0;

while max_error > 0.05 && iterations < max_iterations
    iterations = iterations + 1;
    tic;
    fprintf('FEM Iteration %d with %d total elements\n', iterations, num_elements_per_segment * 3);
    
    % Solve using FEM
    [x_fem, nodal_displacements, element_stresses, error] = solve_fem(A1, A2, A3, E1, E2, L, F1, F2, F3, num_elements_per_segment);
    max_error = max(abs(error));
    
    fprintf('  Maximum error: %.2f%%\n', max_error * 100);
    fprintf('  FEM end displacement: %.4f mm\n', nodal_displacements(end));
    fprintf('  Iteration completed in %.2f seconds\n', toc);
    
    if max_error > 0.05 && iterations < max_iterations
        num_elements_per_segment = num_elements_per_segment * 2;
        fprintf('  Refining mesh to %d elements per segment...\n', num_elements_per_segment);
    end
end

fprintf('\nFinal FEM results:\n');
fprintf('  Number of elements: %d\n', num_elements_per_segment * 3);
fprintf('  Maximum error: %.2f%%\n', max_error * 100);
fprintf('  End displacement: %.4f mm\n', nodal_displacements(end));
fprintf('  Analytical end displacement: %.4f mm\n', total_displacement);
fprintf('  Displacement error: %.2f%%\n\n', abs(nodal_displacements(end) - total_displacement)/total_displacement * 100);

%% 3. VISUALIZATION
fprintf('3. GENERATING VISUALIZATIONS:\n');
fprintf('%s\n', repmat('-', 1, 40));

% Generate enhanced plots in the matlab directory
try
    matlab_dir = 'matlab';
    if ~exist(matlab_dir, 'dir')
        mkdir(matlab_dir);
    end
    enhanced_plotting_matlab(A1, A2, A3, E1, E2, L, F1, F2, F3, num_elements_per_segment, matlab_dir);
    fprintf('Enhanced plots generated and saved to %s successfully.\n', matlab_dir);
catch e
    fprintf('Error generating enhanced plots: %s\n', e.message);
end

%% 4. REPORT GENERATION
fprintf('\n4. GENERATING REPORT:\n');
fprintf('%s\n', repmat('-', 1, 40));

% Generate report
try
    report_file = generate_enhanced_report(A1, A2, A3, E1, E2, L, F1, F2, F3, num_elements_per_segment, matlab_dir);
    fprintf('Report generated: %s\n\n', report_file);
catch e
    fprintf('Error generating report: %s\n', e.message);
end

fprintf('Analysis complete.\n');
