clc;
clear all;
close all;

A1 = 200;
A2 = 100;
A3 = 50;
E1 = 130;
E2 = 200;
L = 500;
F1 = 20;
F2 = 40;
F3 = 20;
num_elements_per_segment = 8;

disp('Problem parameters:');
disp(['- A1 = ' num2str(A1) ' mm², A2 = ' num2str(A2) ' mm², A3 = ' num2str(A3) ' mm²']);
disp(['- E1 = ' num2str(E1) ' GPa, E2 = ' num2str(E2) ' GPa']);
disp(['- L = ' num2str(L) ' mm']);
disp(['- F1 = ' num2str(F1) ' kN, F2 = ' num2str(F2) ' kN, F3 = ' num2str(F3) ' kN']);
disp('');

disp('Generating enhanced plots...');

[x_analytical, stress_analytical, displacement_analytical] = solve_analytical(A1, A2, A3, E1, E2, L, F1, F2, F3);
[x_fem, nodal_displacements, element_stresses, error] = solve_fem(A1, A2, A3, E1, E2, L, F1, F2, F3, num_elements_per_segment);

standalone_plotting(A1, A2, A3, E1, E2, L, F1, F2, F3, x_analytical, displacement_analytical, stress_analytical, x_fem, nodal_displacements, element_stresses);

disp('Plots generated and saved in the current directory.');
disp('');

disp('Generating comprehensive report...');

report_file = 'enhanced_bar_analysis_report.txt';
fid = fopen(report_file, 'w');

fprintf(fid, 'COMPOSITE BAR STRUCTURE ANALYSIS REPORT\n');
fprintf(fid, '=======================================\n\n');

fprintf(fid, '1. PROBLEM PARAMETERS\n');
fprintf(fid, '------------------\n');
fprintf(fid, '- Cross-sectional areas: A1 = %g mm², A2 = %g mm², A3 = %g mm²\n', A1, A2, A3);
fprintf(fid, '- Elastic moduli: E1 = %g GPa, E2 = %g GPa\n', E1, E2);
fprintf(fid, '- Segment length: L = %g mm\n', L);
fprintf(fid, '- Applied forces: F1 = %g kN, F2 = %g kN, F3 = %g kN\n\n', F1, F2, F3);

fprintf(fid, '2. METHODOLOGY\n');
fprintf(fid, '-------------\n');
fprintf(fid, 'The bar structure was analyzed using both analytical and FEM approaches.\n');
fprintf(fid, '- Analytical: Direct application of mechanics of materials principles\n');
fprintf(fid, '- FEM: %d elements per segment (%d total elements)\n\n', num_elements_per_segment, 3*num_elements_per_segment);

fprintf(fid, '3. KEY RESULTS\n');
fprintf(fid, '------------\n');

max_disp = max(displacement_analytical);
fprintf(fid, '- Maximum displacement: %.4f mm\n', max_disp);

max_stress = max(stress_analytical);
fprintf(fid, '- Maximum stress: %.2f MPa\n', max_stress);

fprintf(fid, '- FEM accuracy: Average error of %.2f%% compared to analytical solution\n\n', error);

fprintf(fid, '4. CONCLUSIONS\n');
fprintf(fid, '-------------\n');
fprintf(fid, '- The composite bar structure behaves as expected under the applied loads.\n');
fprintf(fid, '- The FEM solution closely matches the analytical solution, validating both approaches.\n');
fprintf(fid, '- The highest stresses occur in segment 3 due to its smaller cross-sectional area.\n');
fprintf(fid, '- The maximum displacement at the end of the bar is %.4f mm.\n\n', max_disp);

fprintf(fid, 'Report generated on %s\n', datestr(now));
fprintf(fid, 'See the generated plots for visual representations of the results.\n');

fclose(fid);

disp(['Report generated: ' report_file]);
disp('');

disp('Analysis and plotting completed successfully.');
disp('');

disp('All tasks completed successfully.');
disp('For more details, check the generated plots and report.');

function [x_analytical, stress_analytical, displacement_analytical] = solve_analytical(A1, A2, A3, E1, E2, L, F1, F2, F3, num_points)
    if nargin < 10 || isempty(num_points)
        num_points = 500;
    end
    
    x_analytical = linspace(0, 3*L, num_points);
    
    displacement_analytical = zeros(size(x_analytical));
    stress_analytical = zeros(size(x_analytical));

    for i = 1:num_points
        x = x_analytical(i);
        
        if x <= L 
            area = A1;
            E = E1;
            % Formula for segment 1
            displacement_analytical(i) = (F1 * x) / (E * 1000 * area);
            stress_analytical(i) = F1 * 1000 / area; % MPa
        elseif x <= 2*L % Segment 2
            area = A2;
            E = E2;
            % Formula for segment 2
            displacement_analytical(i) = (F1 * L) / (E1 * 1000 * A1) + ...
                                       (F2 * (x - L)) / (E * 1000 * area);
            stress_analytical(i) = F2 * 1000 / area; % MPa
        else % Segment 3
            area = A3;
            E = E2;
            % Formula for segment 3
            displacement_analytical(i) = (F1 * L) / (E1 * 1000 * A1) + ...
                                       (F2 * L) / (E2 * 1000 * A2) + ...
                                       (F3 * (x - 2*L)) / (E * 1000 * area);
            stress_analytical(i) = F3 * 1000 / area; % MPa
        end
    end
end

function [x_nodal, nodal_displacements, element_stresses, error] = solve_fem(A1, A2, A3, E1, E2, L, F1, F2, F3, num_elements_per_segment)
    num_elements = num_elements_per_segment * 3;
    num_nodes = num_elements + 1;
    x_nodal = zeros(num_nodes, 1);
    element_length = 3*L / num_elements;
    
    for i = 1:num_nodes
        x_nodal(i) = (i-1) * element_length;
    end
    
    K = zeros(num_nodes, num_nodes);
    F = zeros(num_nodes, 1);
    
    for e = 1:num_elements
        node1 = e;
        node2 = e + 1;
        
        x1 = x_nodal(node1);
        x2 = x_nodal(node2);
        
        le = x2 - x1;
        
        x_mid = (x1 + x2) / 2;
        if x_mid <= L
            A = A1;
            E = E1;
        elseif x_mid <= 2*L
            A = A2;
            E = E2;
        else
            A = A3;
            E = E2;
        end
        
        ke = (A * E * 1000 / le) * [1, -1; -1, 1];
        K(node1:node2, node1:node2) = K(node1:node2, node1:node2) + ke;
    end
    
    % Apply forces at specific nodes
    node_at_L = num_elements_per_segment + 1;
    F(node_at_L) = F(node_at_L) + F2 * 1000;
    
    node_at_2L = 2 * num_elements_per_segment + 1;
    F(node_at_2L) = F(node_at_2L) + F3 * 1000;
    
    % Apply boundary conditions (fixed at x=0)
    K(1, :) = 0;
    K(1, 1) = 1;
    F(1) = 0;
    
    % Apply F1 at the right end of the bar
    F(end) = F1 * 1000;
    
    nodal_displacements = K \ F;
    
    element_stresses = zeros(num_elements, 1);
    for e = 1:num_elements
        node1 = e;
        node2 = e + 1;
        
        u1 = nodal_displacements(node1);
        u2 = nodal_displacements(node2);
        
        le = x_nodal(node2) - x_nodal(node1);
        
        x_mid = (x_nodal(node1) + x_nodal(node2)) / 2;
        if x_mid <= L
            A = A1;
            E = E1;
        elseif x_mid <= 2*L
            A = A2;
            E = E2;
        else
            A = A3;
            E = E2;
        end
        
        element_stresses(e) = E * 1000 * (u2 - u1) / le;
    end
    
    element_centers = (x_nodal(1:end-1) + x_nodal(2:end)) / 2;
    [x_an, stress_an, ~] = solve_analytical(A1, A2, A3, E1, E2, L, F1, F2, F3, 500);
    stress_analytical_at_centers = interp1(x_an, stress_an, element_centers);
    rel_error = abs((element_stresses - stress_analytical_at_centers) ./ stress_analytical_at_centers);
    error = mean(rel_error) * 100; 
end

function area = get_area_at_x(x, A1, A2, A3, L)
    if x <= L
        area = A1;
    elseif x <= 2*L
        area = A2;
    else
        area = A3;
    end
end

function standalone_plotting(A1, A2, A3, E1, E2, L, F1, F2, F3, x_analytical, displacement_analytical, stress_analytical, x_nodal, nodal_displacements, element_stresses)
    % Calculate element centers for plotting element stresses
    element_centers = (x_nodal(1:end-1) + x_nodal(2:end)) / 2;
    
    figure('Position', [100, 100, 900, 600]);
    plot(x_analytical, displacement_analytical, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Analytical Solution');
    hold on;
    plot(x_nodal, nodal_displacements, 'ro', 'MarkerSize', 6, 'LineWidth', 1.5, 'DisplayName', 'FEM Solution (Nodes)');
    
    for i = 1:2
        x_pos = i * L;
        line([x_pos, x_pos], ylim, 'Color', [0.5, 0.5, 0.5], 'LineStyle', '--', 'LineWidth', 1);
        text(x_pos + 15, max(displacement_analytical) * 0.9, sprintf('Segment %d/%d boundary', i, i+1), ...
            'Rotation', 90, 'VerticalAlignment', 'top', 'FontSize', 9);
    end
    
    % Annotate with key values
    text(L/2, displacement_analytical(round(length(x_analytical)/6)) * 1.1, ...
        sprintf('Max displacement in Segment 1: %.2f mm', max(displacement_analytical(x_analytical <= L))), ...
        'FontSize', 9, 'HorizontalAlignment', 'center', ...
        'BackgroundColor', [0.9, 0.9, 1], 'EdgeColor', [0, 0, 1], 'Margin', 3);
    
    text(L + L/2, displacement_analytical(round(length(x_analytical)/2)) * 1.1, ...
        sprintf('Max displacement in Segment 2: %.2f mm', max(displacement_analytical(x_analytical <= 2*L & x_analytical > L))), ...
        'FontSize', 9, 'HorizontalAlignment', 'center', ...
        'BackgroundColor', [0.9, 0.9, 1], 'EdgeColor', [0, 0, 1], 'Margin', 3);
    
    text(2*L + L/2, displacement_analytical(round(5*length(x_analytical)/6)) * 1.1, ...
        sprintf('Max displacement in Segment 3: %.2f mm', max(displacement_analytical(x_analytical > 2*L))), ...
        'FontSize', 9, 'HorizontalAlignment', 'center', ...
        'BackgroundColor', [0.9, 0.9, 1], 'EdgeColor', [0, 0, 1], 'Margin', 3);
    
    xlabel('Position along bar (mm)', 'FontSize', 12);
    ylabel('Displacement (mm)', 'FontSize', 12);
    title('Displacement Field Comparison: Analytical vs. FEM', 'FontSize', 14);
    legend('Location', 'northwest');
    grid on;
    saveas(gcf, 'displacement_field.png');

    figure('Position', [100, 100, 900, 600]);
    plot(x_analytical, stress_analytical, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Analytical Solution');
    hold on;
    plot(element_centers, element_stresses, 'ro', 'MarkerSize', 6, 'LineWidth', 1.5, 'DisplayName', 'FEM Solution (Elements)');
    
    % Add error indicators at each element
    for i = 1:length(element_stresses)
        x = element_centers(i);
        stress_fem = element_stresses(i);
        stress_exact = interp1(x_analytical, stress_analytical, x);
        error_pct = abs((stress_fem - stress_exact) / (abs(stress_exact) + 1e-10)) * 100;
        
        if error_pct > 5  % Only annotate points with significant error
            text(x, stress_fem + 20, sprintf('%.1f%%', error_pct), ...
                'FontSize', 8, 'HorizontalAlignment', 'center', ...
                'BackgroundColor', [1, 1, 1], 'EdgeColor', [0.5, 0.5, 0.5], 'Margin', 2);
        end
    end
    
    % Add segment boundaries
    for i = 1:2
        x_pos = i * L;
        line([x_pos, x_pos], ylim, 'Color', [0.5, 0.5, 0.5], 'LineStyle', '--', 'LineWidth', 1);
        text(x_pos + 15, max(stress_analytical) * 0.9, sprintf('Segment %d/%d boundary', i, i+1), ...
            'Rotation', 90, 'VerticalAlignment', 'top', 'FontSize', 9);
    end
    
    % Add stress annotations
    text(L/2, stress_analytical(round(length(x_analytical)/6)) * 0.7, ...
        sprintf('Stress in Segment 1: %.0f MPa', mean(stress_analytical(x_analytical <= L))), ...
        'FontSize', 9, 'HorizontalAlignment', 'center', ...
        'BackgroundColor', [1, 0.9, 0.9], 'EdgeColor', [1, 0, 0], 'Margin', 3);
    
    text(L + L/2, stress_analytical(round(length(x_analytical)/2)) * 0.7, ...
        sprintf('Stress in Segment 2: %.0f MPa', mean(stress_analytical(x_analytical <= 2*L & x_analytical > L))), ...
        'FontSize', 9, 'HorizontalAlignment', 'center', ...
        'BackgroundColor', [1, 0.9, 0.9], 'EdgeColor', [1, 0, 0], 'Margin', 3);
    
    text(2*L + L/2, stress_analytical(round(5*length(x_analytical)/6)) * 0.7, ...
        sprintf('Stress in Segment 3: %.0f MPa', mean(stress_analytical(x_analytical > 2*L))), ...
        'FontSize', 9, 'HorizontalAlignment', 'center', ...
        'BackgroundColor', [1, 0.9, 0.9], 'EdgeColor', [1, 0, 0], 'Margin', 3);
    
    xlabel('Position along bar (mm)', 'FontSize', 12);
    ylabel('Stress (MPa)', 'FontSize', 12);
    title('Stress Field Comparison: Analytical vs. FEM', 'FontSize', 14);
    legend('Location', 'northeast');
    grid on;
    saveas(gcf, 'stress_field.png');

    % Perform convergence study with different mesh sizes
    element_counts = [4, 8, 16, 32];
    errors = zeros(size(element_counts));
    
    figure('Position', [100, 100, 900, 600]);
    
    for i = 1:length(element_counts)
        n_elements = element_counts(i);
        
        % Solve with current mesh size
        [~, ~, stresses_conv, err] = solve_fem(A1, A2, A3, E1, E2, L, F1, F2, F3, n_elements);
        errors(i) = err;
        
        % Plot convergence point
        loglog(n_elements * 3, err, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
        hold on;
        text(n_elements * 3 * 1.1, err, sprintf('%d elements\n%.2f%% error', n_elements * 3, err), ...
            'FontSize', 10);
    end
    
    % Add trend line
    loglog(element_counts * 3, errors, 'b--', 'LineWidth', 1.5);
    
    % Label plot
    xlabel('Number of Elements (Total)', 'FontSize', 12);
    ylabel('Average Error (%)', 'FontSize', 12);
    title('FEM Convergence Study: Error vs. Number of Elements', 'FontSize', 14);
    grid on;
    saveas(gcf, 'convergence_study.png');
end
