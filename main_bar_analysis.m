%% Bar Structure Analysis - Main Script
% This script analyzes a composite bar structure using both analytical and
% finite element methods.

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

%% Display parameters
fprintf('Bar Structure Analysis\n');
fprintf('=====================\n');
fprintf('Parameters:\n');
fprintf('A1 = %.0f mm², A2 = %.0f mm², A3 = %.0f mm²\n', A1, A2, A3);
fprintf('E1 = %.1f GPa, E2 = %.1f GPa\n', E1/1000, E2/1000);
fprintf('L = %.0f mm\n', L);
fprintf('F1 = %.1f kN, F2 = %.1f kN, F3 = %.1f kN\n\n', F1/1000, F2/1000, F3/1000);

%% Initial number of elements per segment
num_elements_per_segment = 4;

%% Calculate analytical solution
fprintf('Computing analytical solution...\n');
[x_analytical, stress_analytical, displacement_analytical] = solve_analytical(A1, A2, A3, E1, E2, L, F1, F2, F3);
fprintf('Analytical solution computed.\n\n');

%% Use finite element method with different element counts until error is below 5%
max_error = inf;
max_iterations = 5;
iterations = 0;

while max_error > 0.05 && iterations < max_iterations
    iterations = iterations + 1;
    fprintf('FEM Iteration %d with %d total elements\n', iterations, num_elements_per_segment * 3);
    
    % Solve using FEM
    [x_fem, nodal_displacements, element_stresses, error] = solve_fem(A1, A2, A3, E1, E2, L, F1, F2, F3, num_elements_per_segment);
    max_error = max(abs(error));
    
    fprintf('Maximum error: %.2f%%\n', max_error * 100);
    
    if max_error > 0.05
        num_elements_per_segment = num_elements_per_segment * 2;
    end
end

%% Plot results
fprintf('\nGenerating plots...\n');

% Create plots directory if it doesn't exist
if ~exist('plots', 'dir')
    mkdir('plots');
end

% Element center coordinates for plotting
x_element_centers = (x_fem(1:end-1) + x_fem(2:end)) / 2;

% Figure 1: Displacement Field
figure('Position', [100, 100, 800, 500]);
plot(x_analytical, displacement_analytical, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Analytical');
hold on;
plot(x_fem, nodal_displacements, 'ro-', 'LineWidth', 1.5, 'DisplayName', 'FEM');
title('Displacement Field');
xlabel('Position x (mm)');
ylabel('Displacement (mm)');
grid on;
legend('Location', 'best');
saveas(gcf, 'plots/displacement_field.png');

% Figure 2: Stress Field
figure('Position', [100, 100, 800, 500]);
plot(x_analytical, stress_analytical, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Analytical');
hold on;
plot(x_element_centers, element_stresses, 'ro-', 'LineWidth', 1.5, 'DisplayName', 'FEM');
title('Stress Field');
xlabel('Position x (mm)');
ylabel('Stress (N/mm²)');
grid on;
legend('Location', 'best');
saveas(gcf, 'plots/stress_field.png');

% Figure 3: Error
figure('Position', [100, 100, 800, 500]);
plot(x_element_centers, error, 'bo-', 'LineWidth', 1.5);
hold on;
yline(0.05, 'r--', 'LineWidth', 1.5, 'DisplayName', '5% Error Threshold');
yline(-0.05, 'r--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
title('Relative Error in Stress');
xlabel('Position x (mm)');
ylabel('Relative Error');
grid on;
legend('Location', 'best');
saveas(gcf, 'plots/error.png');

% Figure 4: Area Distribution
figure('Position', [100, 100, 800, 500]);
x_plot = linspace(0, 3*L, 1000);
area = zeros(size(x_plot));
for i = 1:length(x_plot)
    area(i) = get_area_at_x(x_plot(i), A1, A2, A3, L);
end
plot(x_plot, area, 'LineWidth', 1.5);
hold on;
xline(L, 'r--', 'LineWidth', 1.5);
xline(2*L, 'r--', 'LineWidth', 1.5);
title('Cross-sectional Area Distribution');
xlabel('Position x (mm)');
ylabel('Area (mm²)');
grid on;
saveas(gcf, 'plots/area_distribution.png');

fprintf('Plots generated successfully.\n');

%% Generate report
fprintf('Generating report...\n');
generate_report(A1, A2, A3, E1, E2, L, F1, F2, F3, num_elements_per_segment);
fprintf('Report generated: bar_analysis_report.txt\n\n');

fprintf('Analysis complete.\n');
