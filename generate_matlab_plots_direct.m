%% Direct MATLAB Plot Generation Script
% This script directly generates all plots in the matlab directory
% to match the Python implementation exactly

% Clear workspace and figures
clear all;
close all;
clc;

% Output directory
matlab_dir = 'matlab';
if ~exist(matlab_dir, 'dir')
    mkdir(matlab_dir);
end

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

% Define number of elements per segment
num_elements_per_segment = 16;  % Match Python version

fprintf('Generating plots and reports for bar structure analysis...\n');
fprintf('Saving all outputs to %s directory\n\n', matlab_dir);

% Run analytical solution
[x_analytical, stress_analytical, displacement_analytical] = solve_analytical(A1, A2, A3, E1, E2, L, F1, F2, F3);

% Run FEM solution
[x_fem, nodal_displacements, element_stresses, error] = solve_fem(A1, A2, A3, E1, E2, L, F1, F2, F3, num_elements_per_segment);

% Element center coordinates for plotting
x_element_centers = (x_fem(1:end-1) + x_fem(2:end)) / 2;
percent_error = abs(error) * 100;  % Convert to percentage

% Internal force calculations (matches Python implementation)
R = F1 + F2 + F3;  % Total reaction at fixed end
N1 = R;            % Internal force in segment 1
N2 = R - F1;       % Internal force in segment 2
N3 = R - F1 - F2;  % Internal force in segment 3

%% 1. DISPLACEMENT FIELD PLOT
fprintf('Generating enhanced displacement field plot...\n');

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

for i = 1:3
    pos = mid_segment_positions(i);
    prop = segment_properties{i};
    text(pos, max(displacement_analytical)*0.85, prop, ...
         'HorizontalAlignment', 'center', ...
         'BackgroundColor', [1 1 1], 'EdgeColor', [0.7 0.7 0.7], 'Margin', 3);
end

% Add end displacement annotation
text(3*L-100, displacement_analytical(end)*0.7, ...
     sprintf('End displacement\n%.3f mm', displacement_analytical(end)), ...
     'BackgroundColor', [1 1 1], 'EdgeColor', [0.7 0.7 0.7], 'Margin', 3);
% Draw arrow to end point
arrow([3*L-100, displacement_analytical(end)*0.7], [3*L, displacement_analytical(end)]);

% Add grid and styling
grid on;
legend('Location', 'best', 'FontSize', 12);
title('Axial Displacement Field in Composite Bar', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Position along bar, x (mm)', 'FontSize', 12);
ylabel('Axial displacement, u(x) (mm)', 'FontSize', 12);

% Save enhanced displacement plot
saveas(gcf, fullfile(matlab_dir, 'enhanced_displacement_field.png'));

%% 2. STRESS FIELD PLOT
fprintf('Generating enhanced stress field plot...\n');

figure('Position', [100, 100, 900, 700]);
plot(x_analytical, stress_analytical, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Analytical Solution');
hold on;
plot(x_element_centers, element_stresses, 'ro', 'MarkerSize', 6, 'DisplayName', 'FEM Solution (Elements)');

% Add segment boundaries
for i = 1:2
    x_pos = i * L;
    xline(x_pos, '--', 'Color', [0.5 0.5 0.5], 'Alpha', 0.6);
    text(x_pos+15, min(stress_analytical)*1.05, sprintf('Segment %d/%d boundary', i, i+1), ...
         'Rotation', 90, 'VerticalAlignment', 'middle');
end

% Add stress value annotations for each segment
segments = [{0, L, 'σ₁'}, {L, 2*L, 'σ₂'}, {2*L, 3*L, 'σ₃'}];
for s = 1:3
    segment = segments{s};
    start = segment{1};
    ending = segment{2};
    label = segment{3};
    
    % Find indices for this segment
    segment_indices = find(x_analytical >= start & x_analytical <= ending);
    segment_stress = stress_analytical(segment_indices);
    segment_x = x_analytical(segment_indices);
    
    % Add annotation at middle of segment
    mid_idx = floor(length(segment_indices)/2);
    mid_x = segment_x(mid_idx);
    mid_stress = segment_stress(mid_idx);
    
    if start == 0  % First segment has varying stress
        text(mid_x, mid_stress + 100, sprintf('%s: %.1f MPa (varies)', label, mid_stress), ...
             'HorizontalAlignment', 'center', ...
             'BackgroundColor', [1 1 1], 'EdgeColor', [0.7 0.7 0.7], 'Margin', 2);
        arrow([mid_x, mid_stress + 90], [mid_x, mid_stress + 10]);
    else  % Other segments have constant stress
        text(mid_x, mid_stress + 100, sprintf('%s: %.1f MPa (constant)', label, mid_stress), ...
             'HorizontalAlignment', 'center', ...
             'BackgroundColor', [1 1 1], 'EdgeColor', [0.7 0.7 0.7], 'Margin', 2);
        arrow([mid_x, mid_stress + 90], [mid_x, mid_stress + 10]);
    end
end

% Add force arrows
forces = [{L, F1}, {2*L, F2}, {3*L, F3}];
for f = 1:3
    force_data = forces{f};
    pos = force_data{1};
    force = force_data{2};
    
    % Arrow pointing left for force
    arrow([pos, max(stress_analytical)*0.6], [pos-20, max(stress_analytical)*0.6], 'Color', 'r', 'Width', 2);
    
    % Label for force value
    text(pos, max(stress_analytical)*0.7, sprintf('F = %.1f kN', force/1000), ...
         'HorizontalAlignment', 'center', ...
         'BackgroundColor', [1 1 0.8], 'EdgeColor', [1 0.5 0], 'Margin', 3);
end

% Add grid and styling
grid on;
legend('Location', 'best', 'FontSize', 12);
title('Axial Stress Distribution in Composite Bar', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Position along bar, x (mm)', 'FontSize', 12);
ylabel('Axial stress, σ(x) (MPa)', 'FontSize', 12);

% Save enhanced stress plot
saveas(gcf, fullfile(matlab_dir, 'enhanced_stress_field.png'));

%% 3. ERROR DISTRIBUTION PLOT
fprintf('Generating enhanced error plot...\n');

figure('Position', [100, 100, 900, 700]);
bar(x_element_centers, percent_error, 0.7, 'FaceColor', [0.5 0.8 1], 'EdgeColor', [0 0 0.8]);
hold on;
yline(5, '--r', 'LineWidth', 2, 'DisplayName', '5% Error Threshold');

% Add segment boundaries
for i = 1:2
    x_pos = i * L;
    xline(x_pos, '--', 'Color', [0.5 0.5 0.5], 'Alpha', 0.6);
    text(x_pos+15, max(percent_error)*0.5, sprintf('Segment %d/%d boundary', i, i+1), ...
         'Rotation', 90, 'VerticalAlignment', 'middle');
end

% Annotate maximum error
[max_error_value, max_error_idx] = max(percent_error);
max_error_x = x_element_centers(max_error_idx);
text(max_error_x + 100, max_error_value + 1, sprintf('Max error: %.2f%%', max_error_value), ...
     'BackgroundColor', [1 1 1], 'EdgeColor', [0.7 0.7 0.7], 'Margin', 2);
arrow([max_error_x + 90, max_error_value + 1], [max_error_x + 10, max_error_value]);

% Add error statistics
avg_error = mean(percent_error);
std_error = std(percent_error);
text(0.02*3*L, 0.92*max(percent_error), ...
     sprintf('Average error: %.2f%%\nStd deviation: %.2f%%', avg_error, std_error), ...
     'BackgroundColor', [1 1 1], 'EdgeColor', [0.7 0.7 0.7], 'Margin', 3);

% Color bars based on error magnitude
hold on;
for i = 1:length(percent_error)
    if percent_error(i) > 5
        bar(x_element_centers(i), percent_error(i), 0.7, 'FaceColor', 'r', 'EdgeColor', [0.7 0 0]);
    end
end

% Add grid and styling
grid on;
legend('Location', 'best', 'FontSize', 12);
title('Relative Error Between Analytical and FEM Solutions', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Position along bar, x (mm)', 'FontSize', 12);
ylabel('Relative error in stress (%)', 'FontSize', 12);

% Save enhanced error plot
saveas(gcf, fullfile(matlab_dir, 'enhanced_error.png'));

%% 4. AREA DISTRIBUTION PLOT
fprintf('Generating enhanced area distribution plot...\n');

figure('Position', [100, 100, 900, 700]);
x_plot = linspace(0, 3*L, 1000);
area = zeros(size(x_plot));
for i = 1:length(x_plot)
    area(i) = get_area_at_x(x_plot(i), A1, A2, A3, L);
end

plot(x_plot, area, 'g-', 'LineWidth', 3, 'DisplayName', 'Cross-sectional Area');
hold on;

% Add segment boundaries
for i = 1:2
    x_pos = i * L;
    xline(x_pos, '--', 'Color', [0.5 0.5 0.5], 'Alpha', 0.6);
    text(x_pos+15, min(area)*0.9, sprintf('Segment %d/%d boundary', i, i+1), ...
         'Rotation', 90, 'VerticalAlignment', 'middle');
end

% Annotate area values at key positions
text(50, A1 + 10, sprintf('A₁ = %d mm²', A1), 'BackgroundColor', [1 1 1], 'EdgeColor', [0 0.7 0], 'Margin', 2);
arrow([50, A1 + 10], [10, A1]);

text(1.5*L, A2 + 20, sprintf('A₂ = %d mm²', A2), ...
     'HorizontalAlignment', 'center', ...
     'BackgroundColor', [1 1 1], 'EdgeColor', [0 0.7 0], 'Margin', 2);
arrow([1.5*L, A2 + 15], [1.5*L, A2 + 2]);

text(2.5*L, A3 + 20, sprintf('A₃ = %d mm²', A3), ...
     'HorizontalAlignment', 'center', ...
     'BackgroundColor', [1 1 1], 'EdgeColor', [0 0.7 0], 'Margin', 2);
arrow([2.5*L, A3 + 15], [2.5*L, A3 + 2]);

% Add description of the variable section
text(L/2, (A1+A2)/2 + 30, 'Linearly varying\ncross-section', ...
     'HorizontalAlignment', 'center', ...
     'BackgroundColor', [1 1 1], 'EdgeColor', [0 0.7 0], 'Margin', 2);
arrow([L/2, (A1+A2)/2 + 25], [L/2, (A1+A2)/2]);

text(1.5*L - 100, A2 - 40, 'Constant\ncross-section', ...
     'BackgroundColor', [1 1 1], 'EdgeColor', [0 0.7 0], 'Margin', 2);
arrow([1.5*L - 90, A2 - 35], [1.5*L, A2]);

text(2.5*L + 100, A3 - 20, 'Constant\ncross-section', ...
     'BackgroundColor', [1 1 1], 'EdgeColor', [0 0.7 0], 'Margin', 2);
arrow([2.5*L + 90, A3 - 15], [2.5*L, A3]);

% Add material indicators
text(L/2, min(area)*0.7, sprintf('Material 1: E₁ = %.0f GPa', E1/1000), ...
     'HorizontalAlignment', 'center', ...
     'BackgroundColor', [0.8 0.9 1], 'EdgeColor', [0 0 0.7], 'Margin', 2);

text(2.5*L, min(area)*0.7, sprintf('Material 2: E₂ = %.0f GPa', E2/1000), ...
     'HorizontalAlignment', 'center', ...
     'BackgroundColor', [0.8 0.9 1], 'EdgeColor', [0 0 0.7], 'Margin', 2);

% Add grid and styling
grid on;
legend('Location', 'best', 'FontSize', 12);
title('Cross-sectional Area Distribution Along Bar', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Position along bar, x (mm)', 'FontSize', 12);
ylabel('Cross-sectional area, A(x) (mm²)', 'FontSize', 12);

% Save enhanced area distribution plot
saveas(gcf, fullfile(matlab_dir, 'enhanced_area_distribution.png'));

%% 5. COMBINED VISUALIZATION
fprintf('Generating enhanced combined visualization...\n');

figure('Position', [100, 100, 1000, 800]);

% Create a schematic diagram of the bar structure in top subplot
subplot(2, 1, 1);

% Draw the bar baseline
plot([0, 3*L], [0, 0], 'k-', 'LineWidth', 3);
hold on;

% Draw the varying cross-sections
bar_height = 50;

% Segment 1 (varying area)
x_seg1 = linspace(0, L, 100);
width_seg1_top = linspace(A1/10, A2/10, 100);
width_seg1_bottom = -width_seg1_top;
fill(x_seg1, width_seg1_top, 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill(x_seg1, width_seg1_bottom, 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Segment 2 (constant area)
x_seg2 = linspace(L, 2*L, 100);
width_seg2 = A2/10;
fill([L, 2*L, 2*L, L], [width_seg2, width_seg2, -width_seg2, -width_seg2], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Segment 3 (constant area)
x_seg3 = linspace(2*L, 3*L, 100);
width_seg3 = A3/10;
fill([2*L, 3*L, 3*L, 2*L], [width_seg3, width_seg3, -width_seg3, -width_seg3], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Add force arrows
arrow_length = bar_height * 3;
quiver(L, 0, -arrow_length/6, 0, 'r', 'LineWidth', 2, 'MaxHeadSize', 1);
quiver(2*L, 0, -arrow_length/6, 0, 'r', 'LineWidth', 2, 'MaxHeadSize', 1);
quiver(3*L, 0, -arrow_length/6, 0, 'r', 'LineWidth', 2, 'MaxHeadSize', 1);

% Add force labels
text(L-arrow_length/12, 20, sprintf('F₁ = %.0f kN', F1/1000), 'HorizontalAlignment', 'center', ...
     'BackgroundColor', [1 1 1], 'EdgeColor', [1 0 0], 'Margin', 2);
text(2*L-arrow_length/12, 20, sprintf('F₂ = %.0f kN', F2/1000), 'HorizontalAlignment', 'center', ...
     'BackgroundColor', [1 1 1], 'EdgeColor', [1 0 0], 'Margin', 2);
text(3*L-arrow_length/12, 20, sprintf('F₃ = %.0f kN', F3/1000), 'HorizontalAlignment', 'center', ...
     'BackgroundColor', [1 1 1], 'EdgeColor', [1 0 0], 'Margin', 2);

% Add fixed support symbol
plot([-10, 0], [bar_height, bar_height], 'k-', 'LineWidth', 2);
plot([-10, 0], [-bar_height, -bar_height], 'k-', 'LineWidth', 2);
for i = -bar_height:10:bar_height
    plot([-10, 0], [i, 0], 'k-', 'LineWidth', 1);
end

% Add segment labels
text(L/2, -bar_height*2, 'Segment 1\nVarying Area\nE₁', 'HorizontalAlignment', 'center', ...
     'BackgroundColor', [1 1 1], 'EdgeColor', [0.7 0.7 0.7], 'Margin', 2);
text(1.5*L, -bar_height*2, 'Segment 2\nConstant Area\nE₁', 'HorizontalAlignment', 'center', ...
     'BackgroundColor', [1 1 1], 'EdgeColor', [0.7 0.7 0.7], 'Margin', 2);
text(2.5*L, -bar_height*2, 'Segment 3\nConstant Area\nE₂', 'HorizontalAlignment', 'center', ...
     'BackgroundColor', [1 1 1], 'EdgeColor', [0.7 0.7 0.7], 'Margin', 2);

title('Composite Bar Structure - Problem Schematic', 'FontSize', 14, 'FontWeight', 'bold');
axis equal
axis off

% Plot the stress distribution below the schematic
subplot(2, 1, 2);
plot(x_analytical, stress_analytical, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Analytical Stress');
hold on;
plot(x_analytical, displacement_analytical*100, 'g-', 'LineWidth', 2.5, 'DisplayName', 'Analytical Displacement (×100)');
plot(x_element_centers, element_stresses, 'ro', 'MarkerSize', 6, 'DisplayName', 'FEM Stress');
plot(x_fem, nodal_displacements*100, 'mo', 'MarkerSize', 4, 'DisplayName', 'FEM Displacement (×100)');

% Add convergence information
mesh_description = sprintf('%d elements/segment (%d total)', num_elements_per_segment, num_elements_per_segment*3);
error_description = sprintf('Max error: %.2f%%, Avg: %.2f%%', max(abs(error))*100, mean(abs(error))*100);
text(0.98*3*L, 0.05*max(stress_analytical), sprintf('Mesh: %s\n%s', mesh_description, error_description), ...
     'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
     'BackgroundColor', [1 1 1], 'EdgeColor', [0.7 0.7 0.7], 'Margin', 2);

% Add segment boundaries
for i = 1:2
    x_pos = i * L;
    xline(x_pos, '--', 'Color', [0.5 0.5 0.5], 'Alpha', 0.6);
end

% Add grid and styling
grid on;
legend('Location', 'best', 'FontSize', 10);
title('Combined Stress and Displacement Plot (Analytical vs. FEM)', 'FontSize', 14);
xlabel('Position along bar, x (mm)', 'FontSize', 12);
ylabel('Stress (MPa) / Scaled Displacement', 'FontSize', 12);

% Save combined visualization
saveas(gcf, fullfile(matlab_dir, 'enhanced_combined_visualization.png'));

%% 6. Generate report in the matlab directory
fprintf('Generating enhanced report in the matlab directory...\n');

% Generate a basic text report file
report_file = fullfile(matlab_dir, 'matlab_bar_analysis_report.txt');
fid = fopen(report_file, 'w');

fprintf(fid, 'COMPOSITE BAR STRUCTURE ANALYSIS REPORT (MATLAB)\n');
fprintf(fid, '==============================================\n\n');

% Parameters
fprintf(fid, 'PARAMETERS:\n');
fprintf(fid, '  A1 = %.0f mm², A2 = %.0f mm², A3 = %.0f mm²\n', A1, A2, A3);
fprintf(fid, '  E1 = %.1f GPa, E2 = %.1f GPa\n', E1/1000, E2/1000);
fprintf(fid, '  L = %.0f mm\n', L);
fprintf(fid, '  F1 = %.1f kN, F2 = %.1f kN, F3 = %.1f kN\n\n', F1/1000, F2/1000, F3/1000);

% Brief results summary
fprintf(fid, 'RESULTS SUMMARY:\n');
fprintf(fid, '  Total reaction at fixed end: R = %.1f kN\n', R/1000);
fprintf(fid, '  Internal forces: N1 = %.1f kN, N2 = %.1f kN, N3 = %.1f kN\n', N1/1000, N2/1000, N3/1000);
fprintf(fid, '  End displacement: %.4f mm\n', displacement_analytical(end));
fprintf(fid, '  Maximum stress: %.1f MPa\n', max(stress_analytical));
fprintf(fid, '  FEM max error: %.2f%%\n\n', max(percent_error));

% List of generated plots
fprintf(fid, 'GENERATED PLOTS:\n');
fprintf(fid, '  1. enhanced_displacement_field.png\n');
fprintf(fid, '  2. enhanced_stress_field.png\n');
fprintf(fid, '  3. enhanced_error.png\n');
fprintf(fid, '  4. enhanced_area_distribution.png\n');
fprintf(fid, '  5. enhanced_combined_visualization.png\n\n');

fprintf(fid, 'Report generated: %s\n', datestr(now));
fclose(fid);

% Create a markdown file for detailed analysis
md_file = fullfile(matlab_dir, 'matlab_detailed_analysis.md');
fid = fopen(md_file, 'w');

fprintf(fid, '# Bar Structure Analysis - MATLAB Implementation\n\n');
fprintf(fid, 'This analysis was performed using the MATLAB implementation, which matches\n');
fprintf(fid, 'the Python version''s functionality and produces identical visualizations.\n\n');

fprintf(fid, '## Problem Parameters\n\n');
fprintf(fid, '- A1 = %.0f mm², A2 = %.0f mm², A3 = %.0f mm²\n', A1, A2, A3);
fprintf(fid, '- E1 = %.1f GPa, E2 = %.1f GPa\n', E1/1000, E2/1000);
fprintf(fid, '- L = %.0f mm\n', L);
fprintf(fid, '- F1 = %.1f kN, F2 = %.1f kN, F3 = %.1f kN\n\n', F1/1000, F2/1000, F3/1000);

fprintf(fid, '## Generated Visualizations\n\n');
fprintf(fid, 'The MATLAB implementation generates the following visualizations, which match\n');
fprintf(fid, 'the Python implementation exactly:\n\n');
fprintf(fid, '1. Displacement Field (`enhanced_displacement_field.png`)\n');
fprintf(fid, '2. Stress Field (`enhanced_stress_field.png`)\n');
fprintf(fid, '3. Error Distribution (`enhanced_error.png`)\n');
fprintf(fid, '4. Area Distribution (`enhanced_area_distribution.png`)\n');
fprintf(fid, '5. Combined Visualization (`enhanced_combined_visualization.png`)\n\n');

fprintf(fid, 'All visualizations include detailed annotations, force arrows, and proper formatting\n');
fprintf(fid, 'to match the Python implementation.\n\n');

fclose(fid);

fprintf('\nAll plots and reports generated successfully in the %s directory.\n', matlab_dir);
fprintf('The following files were created:\n');
fprintf('  - enhanced_displacement_field.png\n');
fprintf('  - enhanced_stress_field.png\n');
fprintf('  - enhanced_error.png\n');
fprintf('  - enhanced_area_distribution.png\n');
fprintf('  - enhanced_combined_visualization.png\n');
fprintf('  - matlab_bar_analysis_report.txt\n');
fprintf('  - matlab_detailed_analysis.md\n');
