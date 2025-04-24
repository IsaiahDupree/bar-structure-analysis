function enhanced_plotting_matlab(A1, A2, A3, E1, E2, L, F1, F2, F3, num_elements_per_segment, save_dir)
% ENHANCED_PLOTTING_MATLAB - Creates improved visualizations with detailed annotations
% This function creates enhanced plots for bar structure analysis with detailed annotations
% and better styling using the refined FEM methodology, exactly matching the Python version
%
% Syntax:
%   enhanced_plotting_matlab(A1, A2, A3, E1, E2, L, F1, F2, F3, num_elements_per_segment, save_dir)
%
% Inputs:
%   A1, A2, A3 - Cross-sectional areas in mm²
%   E1, E2 - Elastic moduli in N/mm²
%   L - Length of each segment in mm
%   F1, F2, F3 - Forces applied at segment ends in N
%   num_elements_per_segment - Number of elements per segment for FEM
%   save_dir - Directory to save the plots (default: 'matlab')
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

%% Handle default save directory
if nargin < 11
    save_dir = 'matlab';
end

% Create a directory for saving plots if it doesn't exist
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

%% Calculate analytical solution
[x_analytical, stress_analytical, displacement_analytical] = solve_analytical(A1, A2, A3, E1, E2, L, F1, F2, F3);

%% Current mesh settings
num_elements_per_segment = 8;  % Base calculation with medium-density mesh

% Solve using FEM
[x_fem, nodal_displacements, element_stresses, error] = solve_fem(A1, A2, A3, E1, E2, L, F1, F2, F3, num_elements_per_segment);

% Element center coordinates for plotting
x_element_centers = (x_fem(1:end-1) + x_fem(2:end)) / 2;
percent_error = abs(error) * 100;  % Convert to percentage

% Calculate additional analytical values for correctness validation
% Similar to the Python implementation's approach
R = F1 + F2 + F3;  % Total reaction at fixed end
N1 = R;            % Internal force in segment 1
N2 = R - F1;       % Internal force in segment 2
N3 = R - F1 - F2;  % Internal force in segment 3

% 1. ENHANCED DISPLACEMENT FIELD PLOT
figure('Position', [100, 100, 900, 700]);
plot(x_analytical, displacement_analytical, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Analytical Solution');
hold on;
plot(x_fem, nodal_displacements, 'ro-', 'MarkerSize', 6, 'LineWidth', 1.5, 'DisplayName', 'FEM Solution (Nodes)');

% Add segment boundaries
for i = 1:2
    x_pos = i * L;
    xline(x_pos, '--', 'Color', [0.5 0.5 0.5], 'Alpha', 0.6);
    text(x_pos+15, max(displacement_analytical)*0.15, sprintf('Segment %d/%d boundary', i, i+1), ...
         'Rotation', 90, 'VerticalAlignment', 'middle');
end

% Add segment annotations
mid_segment_positions = [L/2, L*1.5, L*2.5];
segment_properties = {...
    sprintf('A₁=%d-%d mm², E₁=%.0f GPa', A1, A2, E1/1000), ...
    sprintf('A₂=%d mm², E₁=%.0f GPa', A2, E1/1000), ...
    sprintf('A₃=%d mm², E₂=%.0f GPa', A3, E2/1000) ...
};
xlabel('Position (mm)', 'FontSize', 12);
ylabel('Displacement (mm)', 'FontSize', 12);
title('Displacement Field: Analytical vs. FEM Solutions', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 11);

% Add convergence indicator
annotation('textbox', [0.5, 0.02, 0, 0], 'String', ...
    sprintf('Maximum Error: %.2f%% (with %d elements per segment)', max(percent_error), num_elements_per_segment), ...
    'HorizontalAlignment', 'center', 'FitBoxToText', 'on', 'BackgroundColor', 'white', ...
    'EdgeColor', [0.8, 0.8, 0.8], 'LineWidth', 1, 'FontSize', 10);

% Save the enhanced displacement field plot
saveas(gcf, fullfile(save_dir, 'enhanced_displacement_field.png'));

%% STRESS FIELD PLOT
figure('Position', [100, 100, 900, 600]);
plot(x_analytical, stress_analytical, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Analytical Solution');
hold on;
plot(x_element_centers, element_stresses, 'ro', 'MarkerSize', 6, 'LineWidth', 1.5, 'DisplayName', 'FEM Solution (Elements)');

% Add error indicators at each element
for i = 1:length(element_stresses)
    x = x_element_centers(i);
    stress_fem = element_stresses(i);
    stress_exact = interp1(x_analytical, stress_analytical, x);
    error_pct = abs((stress_fem - stress_exact) / (abs(stress_exact) + 1e-10)) * 100;
    
    if error_pct > 5  % Only annotate points with significant error
        text(x, stress_fem + 20, sprintf('%.1f%%', error_pct), ...
            'FontSize', 8, 'HorizontalAlignment', 'center', ...
            'BackgroundColor', 'white', 'EdgeColor', 'gray', 'Margin', 2);
    end
end

% Add segment boundaries
for i = 1:2
    x_pos = i * L;
    line([x_pos, x_pos], ylim, 'Color', [0.5, 0.5, 0.5], 'LineStyle', '--', 'LineWidth', 1);
    text(x_pos + 15, max(stress_analytical) * 0.9, sprintf('Segment %d/%d boundary', i, i+1), ...
        'Rotation', 90, 'VerticalAlignment', 'top', 'FontSize', 9);
end

% Add material property annotations
annotation('textbox', [0.2, 0.85, 0.1, 0.05], 'String', ...
    sprintf('E₁ = %.0f GPa', E1/1000), 'FitBoxToText', 'on', ...
    'BackgroundColor', 'lightblue', 'EdgeColor', 'gray', 'LineWidth', 1);
annotation('textbox', [0.45, 0.85, 0.1, 0.05], 'String', ...
    sprintf('E₁ = %.0f GPa', E1/1000), 'FitBoxToText', 'on', ...
    'BackgroundColor', 'lightblue', 'EdgeColor', 'gray', 'LineWidth', 1);
annotation('textbox', [0.75, 0.85, 0.1, 0.05], 'String', ...
    sprintf('E₂ = %.0f GPa', E2/1000), 'FitBoxToText', 'on', ...
    'BackgroundColor', 'lightblue', 'EdgeColor', 'gray', 'LineWidth', 1);

% Add grid and styling
grid on;
xlabel('Position (mm)', 'FontSize', 12);
ylabel('Stress (MPa)', 'FontSize', 12);
title('Stress Field: Analytical vs. FEM Solutions', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 11);

% Add accuracy indicator
annotation('textbox', [0.5, 0.02, 0, 0], 'String', ...
    sprintf('Average Error: %.2f%% (with %d elements per segment)', mean(percent_error), num_elements_per_segment), ...
    'HorizontalAlignment', 'center', 'FitBoxToText', 'on', 'BackgroundColor', 'white', ...
    'EdgeColor', [0.8, 0.8, 0.8], 'LineWidth', 1, 'FontSize', 10);

% Save the enhanced stress field plot
saveas(gcf, fullfile(save_dir, 'enhanced_stress_field.png'));

%% CONVERGENCE STUDY PLOT
figure('Position', [100, 100, 800, 500]);
loglog(mesh_sizes, convergence_errors, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
hold on;

% Add theoretical convergence line for reference (slope -1)
ref_x = [mesh_sizes(1), mesh_sizes(end)];
ref_y = convergence_errors(1) * (ref_x(1) ./ ref_x).^1;  % Order of convergence = 1
loglog(ref_x, ref_y, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Theoretical 1st Order Convergence');

% Add data points with values
for i = 1:length(mesh_sizes)
    text(mesh_sizes(i), convergence_errors(i) * 1.1, sprintf('%.2f%%', convergence_errors(i)), ...
        'FontSize', 9, 'HorizontalAlignment', 'center', ...
        'BackgroundColor', 'white', 'EdgeColor', 'gray', 'Margin', 2);
end

% Add grid and styling
grid on;
xlabel('Number of Elements per Segment', 'FontSize', 12);
ylabel('Maximum Error (%)', 'FontSize', 12);
title('Convergence Study: Error vs. Mesh Density', 'FontSize', 14);
legend('Error', 'Theoretical 1st Order Convergence', 'Location', 'best', 'FontSize', 11);

% Save the convergence study plot
saveas(gcf, fullfile(save_dir, 'enhanced_convergence_study.png'));

%% AREA DISTRIBUTION PLOT
figure('Position', [100, 100, 800, 400]);

% Create x points for area calculation
x_area = linspace(0, 3*L, 1000);

% Calculate area at each point
area_values = zeros(size(x_area));
for i = 1:length(x_area)
    x = x_area(i);
    if x <= L
        % First segment: Linear interpolation from A1 to A2
        area_values(i) = A1 + (A2 - A1) * (x / L);
    elseif x <= 2*L
        % Second segment: Constant A2
        area_values(i) = A2;
    else
        % Third segment: Constant A3
        area_values(i) = A3;
    end
end

% Plot area distribution
plot(x_area, area_values, 'k-', 'LineWidth', 2.5);

% Add segment boundaries
for i = 1:2
    x_pos = i * L;
    line([x_pos, x_pos], ylim, 'Color', [0.5, 0.5, 0.5], 'LineStyle', '--', 'LineWidth', 1);
end

% Add grid and styling
grid on;
xlabel('Position (mm)', 'FontSize', 12);
ylabel('Cross-sectional Area (mm²)', 'FontSize', 12);
title('Cross-sectional Area Distribution', 'FontSize', 14);

% Save the area distribution plot
saveas(gcf, fullfile(save_dir, 'enhanced_area_distribution.png'));

%% COMBINED VISUALIZATION WITH BAR SCHEMATIC
figure('Position', [100, 100, 1000, 800]);

% Create a 2-panel plot
subplot(2, 1, 1);

% Draw the bar schematic
% First segment (tapered)
x1 = [0, L];
y_top1 = [A1/20, A2/20];  % Scale down for visualization
y_bottom1 = [-A1/20, -A2/20];
fill([x1, fliplr(x1)], [y_top1, fliplr(y_bottom1)], [0.8, 0.8, 1.0], 'EdgeColor', 'blue', 'LineWidth', 1.5);

% Second segment (constant area)
x2 = [L, 2*L];
y_top2 = [A2/20, A2/20];
y_bottom2 = [-A2/20, -A2/20];
hold on;
fill([x2, fliplr(x2)], [y_top2, fliplr(y_bottom2)], [0.8, 1.0, 0.8], 'EdgeColor', 'green', 'LineWidth', 1.5);

% Third segment (constant area)
x3 = [2*L, 3*L];
y_top3 = [A3/20, A3/20];
y_bottom3 = [-A3/20, -A3/20];
fill([x3, fliplr(x3)], [y_top3, fliplr(y_bottom3)], [1.0, 0.8, 0.8], 'EdgeColor', 'red', 'LineWidth', 1.5);

% Add material labels
text(L/2, A1/15, sprintf('Segment 1\nA₁ = %d mm²\nE₁ = %d GPa', A1, E1/1000), 'HorizontalAlignment', 'center');
text(1.5*L, A2/15, sprintf('Segment 2\nA₂ = %d mm²\nE₁ = %d GPa', A2, E1/1000), 'HorizontalAlignment', 'center');
text(2.5*L, A3/15, sprintf('Segment 3\nA₃ = %d mm²\nE₂ = %d GPa', A3, E2/1000), 'HorizontalAlignment', 'center');

% Add force arrows
arrow_length = 30;
% F1 arrow at x = L
arrow([L, 0], [L+arrow_length, 0], 'Length', 12, 'Width', 5, 'EdgeColor', 'k', 'FaceColor', 'r');
text(L+arrow_length/2, -10, sprintf('F₁ = %d kN', F1/1000), 'HorizontalAlignment', 'center');

% F2 arrow at x = 2L
arrow([2*L, 0], [2*L+arrow_length, 0], 'Length', 12, 'Width', 5, 'EdgeColor', 'k', 'FaceColor', 'r');
text(2*L+arrow_length/2, -10, sprintf('F₂ = %d kN', F2/1000), 'HorizontalAlignment', 'center');

% F3 arrow at x = 3L
arrow([3*L, 0], [3*L+arrow_length, 0], 'Length', 12, 'Width', 5, 'EdgeColor', 'k', 'FaceColor', 'r');
text(3*L+arrow_length/2, -10, sprintf('F₃ = %d kN', F3/1000), 'HorizontalAlignment', 'center');

% Fixed end indication
line([0, 0], [-20, 20], 'Color', 'k', 'LineWidth', 2);
for i = 1:5
    line([0, -10], [-20+8*i, -20+8*i], 'Color', 'k', 'LineWidth', 2);
end
text(-20, 0, 'Fixed end', 'Rotation', 90, 'HorizontalAlignment', 'center');

% Set axis properties
axis equal;
axis([-(arrow_length+20), 3*L+(arrow_length+20), -A1/10, A1/10]);
title('Bar Structure Schematic (Not to Scale)');
xlabel('Position (mm)');
set(gca, 'YTick', []);  % Remove y-axis ticks for schematic

% Plot the stress distribution below the schematic
subplot(2, 1, 2);
plot(x_analytical, stress_analytical, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Analytical Stress');
hold on;
plot(x_analytical, displacement_analytical*100, 'g-', 'LineWidth', 2.5, 'DisplayName', 'Analytical Displacement (×100)');
plot(x_element_centers, element_stresses, 'ro', 'MarkerSize', 6, 'DisplayName', 'FEM Stress');
plot(x_fem, nodal_displacements*100, 'mo', 'MarkerSize', 4, 'DisplayName', 'FEM Displacement (×100)');

% Add convergence information
mesh_description = sprintf('%d elements/segment (%d total)', num_elements_per_segment, num_elements_per_segment*3);
error_description = sprintf('Max error: %.2f%%, Avg: %.2f%%', max(percent_error), mean(percent_error));
annotation('textbox', [0.8, 0.25, 0.15, 0.05], 'String', ...
    {mesh_description, error_description}, ...
    'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'gray', ...
    'LineWidth', 1, 'FontSize', 9);

% Add segment boundaries
for i = 1:2
    x_pos = i * L;
    line([x_pos, x_pos], ylim, 'Color', [0.5, 0.5, 0.5], 'LineStyle', '--', 'LineWidth', 1);
end

% Add grid and styling
grid on;
xlabel('Position (mm)', 'FontSize', 12);
ylabel('Stress (MPa) / Displacement×100 (mm)', 'FontSize', 12);
title('Stress and Displacement Distribution', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 10);

% Save the combined visualization
saveas(gcf, fullfile(save_dir, 'enhanced_combined_visualization.png'));

fprintf('Enhanced plots created and saved to %s directory\n', fullfile(pwd, save_dir));
