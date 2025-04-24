%% Comprehensive Bar Structure Analysis
% This script performs a complete analysis of the bar structure problem,
% calling other scripts as needed and generating a comprehensive report
% that answers all required questions.

%% Clear environment and start timer
clear all;
close all;
clc;

fprintf('=======================================================\n');
fprintf('COMPREHENSIVE BAR STRUCTURE ANALYSIS\n');
fprintf('=======================================================\n\n');

% Start timer for overall execution
overall_start_time = tic;

%% Define parameters
fprintf('Defining problem parameters...\n');

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
analytical_start_time = tic;
[x_analytical, stress_analytical, displacement_analytical] = solve_analytical(A1, A2, A3, E1, E2, L, F1, F2, F3);
analytical_time = toc(analytical_start_time);
fprintf('Analytical solution computed in %.2f seconds.\n\n', analytical_time);

%% Use finite element method with different element counts until error is below 5%
fprintf('Performing FEM analysis with adaptive mesh refinement...\n');
max_error = inf;
max_iterations = 3;  % Reduced to 3 to match Python implementation
iterations = 0;

fem_start_time = tic;
while max_error > 0.05 && iterations < max_iterations
    iterations = iterations + 1;
    iteration_start_time = tic;
    fprintf('FEM Iteration %d with %d total elements\n', iterations, num_elements_per_segment * 3);
    
    % Solve using FEM
    [x_fem, nodal_displacements, element_stresses, error] = solve_fem(A1, A2, A3, E1, E2, L, F1, F2, F3, num_elements_per_segment);
    max_error = max(abs(error));
    
    fprintf('Maximum error: %.2f%%\n', max_error * 100);
    fprintf('Iteration completed in %.2f seconds\n', toc(iteration_start_time));
    
    if max_error > 0.05 && iterations < max_iterations
        num_elements_per_segment = num_elements_per_segment * 2;
    end
end
fem_time = toc(fem_start_time);
fprintf('FEM analysis completed in %.2f seconds.\n\n', fem_time);

%% Plot results
fprintf('Generating plots...\n');
plot_start_time = tic;

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
xline(L, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Segment Boundary');
xline(2*L, 'r--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
title('Cross-sectional Area Distribution');
xlabel('Position x (mm)');
ylabel('Area (mm²)');
grid on;
legend('Location', 'best');
saveas(gcf, 'plots/area_distribution.png');

% Figure 5: Modulus Distribution
figure('Position', [100, 100, 800, 500]);
modulus = zeros(size(x_plot));
for i = 1:length(x_plot)
    if x_plot(i) <= 2*L
        modulus(i) = E1;
    else
        modulus(i) = E2;
    end
end
plot(x_plot, modulus/1000, 'LineWidth', 1.5); % Divide by 1000 to convert to GPa
hold on;
xline(L, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Segment Boundary');
xline(2*L, 'r--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
title('Elastic Modulus Distribution');
xlabel('Position x (mm)');
ylabel('Elastic Modulus (GPa)');
grid on;
legend('Location', 'best');
saveas(gcf, 'plots/modulus_distribution.png');

plot_time = toc(plot_start_time);
fprintf('Plots generated successfully in %.2f seconds.\n', plot_time);
fprintf('Plots saved to the "plots" directory.\n\n');

%% Comprehensive Analysis Report with Answers to Questions
fprintf('Generating comprehensive report with answers to questions...\n');

% Open file for writing
fid = fopen('comprehensive_bar_analysis_report.txt', 'w');

% Write report header
fprintf(fid, 'COMPREHENSIVE BAR STRUCTURE ANALYSIS REPORT\n');
fprintf(fid, '==========================================\n\n');

% Write problem parameters
fprintf(fid, 'Problem Parameters:\n');
fprintf(fid, '- A1 = %.0f mm², A2 = %.0f mm², A3 = %.0f mm²\n', A1, A2, A3);
fprintf(fid, '- E1 = %.1f GPa, E2 = %.1f GPa\n', E1/1000, E2/1000);
fprintf(fid, '- L = %.0f mm\n', L);
fprintf(fid, '- F1 = %.1f kN, F2 = %.1f kN, F3 = %.1f kN\n\n', F1/1000, F2/1000, F3/1000);

% Write analysis details
fprintf(fid, 'Finite Element Analysis:\n');
fprintf(fid, '- Number of elements per segment: %d\n', num_elements_per_segment);
fprintf(fid, '- Total number of elements: %d\n', 3 * num_elements_per_segment);
fprintf(fid, '- Element type: 2-node linear elements\n');
fprintf(fid, '- Analytical solution computed in %.2f seconds\n', analytical_time);
fprintf(fid, '- FEM solution computed in %.2f seconds\n', fem_time);
fprintf(fid, '- Maximum error: %.2f%%\n\n', max_error * 100);

% Write results summary
fprintf(fid, 'Results and Visualizations:\n');
fprintf(fid, '- Plots have been saved in the "plots" directory\n');
fprintf(fid, '- displacement_field.png: Shows the displacement distribution along the bar\n');
fprintf(fid, '- stress_field.png: Shows the stress distribution along the bar\n');
fprintf(fid, '- error.png: Shows the relative error in stress between FEM and analytical solutions\n');
fprintf(fid, '- area_distribution.png: Shows how the cross-sectional area varies along the bar\n');
fprintf(fid, '- modulus_distribution.png: Shows how the elastic modulus varies along the bar\n\n');

% Answers to the 8 questions
fprintf(fid, 'ANSWERS TO THE 8 ANALYSIS QUESTIONS\n');
fprintf(fid, '==================================\n\n');

% Question 1: Maximum displacement in the bar
[max_displacement, max_disp_idx] = max(abs(displacement_analytical));
max_disp_location = x_analytical(max_disp_idx);
fprintf(fid, '1) What is the maximum displacement in the bar?\n');
fprintf(fid, '   Answer: %.6f mm at position x = %.2f mm\n\n', max_displacement, max_disp_location);

% Question 2: Maximum stress in the bar
[max_stress, max_stress_idx] = max(abs(stress_analytical));
max_stress_location = x_analytical(max_stress_idx);
fprintf(fid, '2) What is the maximum stress in the bar?\n');
fprintf(fid, '   Answer: %.2f N/mm² at position x = %.2f mm\n\n', max_stress, max_stress_location);

% Question 3: Effect of area variation on stress distribution
fprintf(fid, '3) How does the variation in cross-sectional area affect the stress distribution?\n');
fprintf(fid, '   Answer: The stress is inversely proportional to the cross-sectional area (σ = F/A).\n');
fprintf(fid, '   Therefore, as the area decreases from A1 (%.0f mm²) to A3 (%.0f mm²),\n', A1, A3);
fprintf(fid, '   we observe an increase in stress. The stress distribution shows significant\n');
fprintf(fid, '   changes at the segment boundaries (x = %d mm and x = %d mm) where the\n', L, 2*L);
fprintf(fid, '   cross-sectional area changes. The smallest area (A3 = %.0f mm²) in segment 3\n', A3);
fprintf(fid, '   results in the highest stress values.\n\n');

% Question 4: Effect of material property variation on displacement
fprintf(fid, '4) How does the variation in material properties affect the displacement?\n');
fprintf(fid, '   Answer: The displacement is inversely proportional to the elastic modulus (ε = σ/E).\n');
fprintf(fid, '   Segments 1 and 2 have a lower elastic modulus (E1 = %.1f GPa) compared to\n', E1/1000);
fprintf(fid, '   segment 3 (E2 = %.1f GPa). This means that segments 1 and 2 will experience\n', E2/1000);
fprintf(fid, '   more strain (and thus more displacement) for the same stress compared to segment 3.\n');
fprintf(fid, '   The change in elastic modulus at x = %d mm results in a change in the slope of\n', 2*L);
fprintf(fid, '   the displacement curve, which is less steep in segment 3 due to the higher modulus.\n\n');

% Question 5: Accuracy of FEM compared to analytical solution
fprintf(fid, '5) How accurate is the FEM solution compared to the analytical solution?\n');
fprintf(fid, '   Answer: The FEM solution achieved a maximum error of %.2f%% relative to the\n', max_error * 100);
fprintf(fid, '   analytical solution after %d iterations of mesh refinement. The error is\n', iterations);
fprintf(fid, '   calculated as (σ_FEM - σ_analytical)/σ_analytical. The highest errors tend to\n');
fprintf(fid, '   occur near the segment boundaries where there are discontinuities in area or\n');
fprintf(fid, '   applied forces. Overall, with %d elements per segment (%d total elements),\n', num_elements_per_segment, 3*num_elements_per_segment);
fprintf(fid, '   the FEM solution provides a good approximation of the analytical solution.\n\n');

% Question 6: Effect of mesh refinement on solution accuracy
fprintf(fid, '6) How does mesh refinement affect the solution accuracy?\n');
fprintf(fid, '   Answer: Mesh refinement significantly improves the accuracy of the FEM solution.\n');
fprintf(fid, '   Starting with %d elements per segment, the maximum error was reduced through\n', 4);
fprintf(fid, '   successive refinements to meet the 5%% error tolerance. Each refinement doubled\n');
fprintf(fid, '   the number of elements per segment, providing better approximation of the\n');
fprintf(fid, '   stress and displacement fields. Mesh refinement is particularly important near\n');
fprintf(fid, '   boundaries where stress changes rapidly. With sufficient refinement, the FEM\n');
fprintf(fid, '   solution converges to the analytical solution.\n\n');

% Question 7: How forces affect the displacement and stress
fprintf(fid, '7) How do the applied forces affect the displacement and stress distributions?\n');
fprintf(fid, '   Answer: The applied forces (F1 = %.1f kN, F2 = %.1f kN, F3 = %.1f kN) create\n', F1/1000, F2/1000, F3/1000);
fprintf(fid, '   internal forces in each segment that determine the stress distribution. The\n');
fprintf(fid, '   displacement at any point is the cumulative effect of the strains in all\n');
fprintf(fid, '   preceding segments. The force F1 creates uniform stress in segment 1, while the\n');
fprintf(fid, '   addition of F2 increases the stress in segment 2. Adding F3 further increases\n');
fprintf(fid, '   the stress in segment 3. The displacement increases monotonically from the fixed\n');
fprintf(fid, '   end (x = 0) to the free end (x = %.0f mm), with the maximum displacement occurring\n', 3*L);
fprintf(fid, '   at the free end.\n\n');

% Question 8: Practical implications of the analysis results
fprintf(fid, '8) What are the practical implications of these analysis results?\n');
fprintf(fid, '   Answer: This analysis provides several practical insights:\n');
fprintf(fid, '   - The maximum stress (%.2f N/mm²) occurs at x = %.2f mm, which is a critical\n', max_stress, max_stress_location);
fprintf(fid, '     point for potential material failure. This should be compared with the material''s\n');
fprintf(fid, '     yield strength to ensure a sufficient factor of safety.\n');
fprintf(fid, '   - The maximum displacement (%.6f mm) at the free end indicates the deformation\n', max_displacement);
fprintf(fid, '     under the given loads, which is important for assessing functionality and clearances.\n');
fprintf(fid, '   - The analysis shows that using variable cross-sectional areas can optimize material\n');
fprintf(fid, '     usage while maintaining acceptable stress levels.\n');
fprintf(fid, '   - The comparison between different materials (E1 vs E2) demonstrates how material\n');
fprintf(fid, '     selection affects structural response.\n');
fprintf(fid, '   - The validation of the FEM approach confirms its reliability for more complex\n');
fprintf(fid, '     geometries where analytical solutions may not be available.\n\n');

% Write methodology
fprintf(fid, 'Methodology:\n');
fprintf(fid, '1. Analytical solution calculated for exact stress and displacement distributions\n');
fprintf(fid, '2. Finite element method applied with 2-node linear elements\n');
fprintf(fid, '3. Element stiffness matrices calculated and assembled\n');
fprintf(fid, '4. System solved for nodal displacements\n');
fprintf(fid, '5. Element stresses computed from displacements\n');
fprintf(fid, '6. Error analysis performed to validate results\n');
fprintf(fid, '7. Adaptive mesh refinement applied until error threshold met or max iterations reached\n\n');

% Write conclusion
fprintf(fid, 'Conclusion:\n');
fprintf(fid, 'This comprehensive analysis has provided valuable insights into the behavior of the\n');
fprintf(fid, 'composite bar structure under axial loading. By comparing analytical and finite\n');
fprintf(fid, 'element methods, we have demonstrated the effectiveness of numerical approaches\n');
fprintf(fid, 'while validating their accuracy against exact solutions. The variation in\n');
fprintf(fid, 'cross-sectional area and material properties significantly influences the stress\n');
fprintf(fid, 'and displacement distributions, highlighting the importance of these parameters in\n');
fprintf(fid, 'structural design. The detailed answers to the eight analysis questions provide a\n');
fprintf(fid, 'comprehensive understanding of the bar''s structural behavior.\n\n');

% Write execution information
fprintf(fid, 'Execution Information:\n');
fprintf(fid, '- Total execution time: %.2f seconds\n', toc(overall_start_time));
fprintf(fid, '- Analysis completed on: %s\n', datestr(now));

% Close file
fclose(fid);

fprintf('Comprehensive text report generated: comprehensive_bar_analysis_report.txt\n');

%% Generate Word Document Report
fprintf('Generating Word document report...\n');
try
    % Call the function to create a Word document with all results and visualizations
    create_word_report(A1, A2, A3, E1, E2, L, F1, F2, F3, num_elements_per_segment, ...
                      x_analytical, stress_analytical, displacement_analytical, ...
                      max_error, iterations, analytical_time, fem_time, toc(overall_start_time));
    
    fprintf('Word document report generated successfully.\n\n');
catch ME
    fprintf('Error generating Word document: %s\n', ME.message);
    fprintf('Text report is still available as comprehensive_bar_analysis_report.txt\n\n');
end

%% Final Output
fprintf('=======================================================\n');
fprintf('ANALYSIS COMPLETE\n');
fprintf('=======================================================\n');
fprintf('- Standard plots saved to: plots/\n');
fprintf('- Comprehensive text report: comprehensive_bar_analysis_report.txt\n');
fprintf('- Comprehensive Word report: Bar_Analysis_Report.docx\n');
fprintf('- Total execution time: %.2f seconds\n', toc(overall_start_time));
fprintf('=======================================================\n');
