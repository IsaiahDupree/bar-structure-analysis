%% Bar Structure Analysis - Complete MATLAB Solution
% This script runs the complete bar structure analysis with enhanced plotting
% All required functions are included in this file

clc;
clear all;
close all;

% Add parent directory to path to find all required functions
parent_dir = fullfile(pwd, '..');
addpath(parent_dir);
disp(['Added parent directory to path: ' parent_dir]);

disp('Bar Structure Analysis - Complete MATLAB Solution');
disp('===============================================');

%% Define parameters for the analysis
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
num_elements_per_segment = 8; % Number of elements per segment

disp('Problem parameters:');
disp(['- A1 = ' num2str(A1) ' mm², A2 = ' num2str(A2) ' mm², A3 = ' num2str(A3) ' mm²']);
disp(['- E1 = ' num2str(E1) ' GPa, E2 = ' num2str(E2) ' GPa']);
disp(['- L = ' num2str(L) ' mm']);
disp(['- F1 = ' num2str(F1) ' kN, F2 = ' num2str(F2) ' kN, F3 = ' num2str(F3) ' kN']);
disp('');

%% Step 1: Generate enhanced plots
disp('Generating enhanced plots...');

% Create plots directory if it doesn't exist
plots_dir = '../plots';
if ~exist(plots_dir, 'dir')
    mkdir(plots_dir);
end

% Run the enhanced plotting function with parameters
enhanced_plotting_matlab(A1, A2, A3, E1, E2, L, F1, F2, F3, num_elements_per_segment, plots_dir);

disp(['Plots saved to ' fullfile(pwd, '../plots')]);
disp('');

disp('All tasks completed successfully.');
disp('For more details, check the generated plots.');

%% Helper Functions

% Main enhanced plotting function
function enhanced_plotting_matlab(A1, A2, A3, E1, E2, L, F1, F2, F3, num_elements_per_segment, plots_dir)
    %% Create shared variables for plotting
    total_length = 3 * L;
    num_points = 500;
    x_analytical = linspace(0, total_length, num_points);
    
    % Calculate analytical solution
    [displacement_analytical, stress_analytical] = analytical_solution(A1, A2, A3, E1, E2, L, F1, F2, F3, x_analytical);
    
    % Calculate FEM solution
    [elements, nodes, nodal_displacements, element_stresses] = ...
        solve_fem(A1, A2, A3, E1, E2, L, F1, F2, F3, num_elements_per_segment);
    
    % Extract data for plotting
    x_fem = nodes;
    u_fem = nodal_displacements;
    
    % Calculate element center positions
    x_element_centers = zeros(length(elements), 1);
    for i = 1:length(elements)
        node1 = elements(i, 1);
        node2 = elements(i, 2);
        x_element_centers(i) = (nodes(node1) + nodes(node2)) / 2;
    end
    
    %% DISPLACEMENT FIELD PLOT
    figure('Position', [100, 100, 900, 600]);
    plot(x_analytical, displacement_analytical, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Analytical Solution');
    hold on;
    plot(x_fem, u_fem, 'ro', 'MarkerSize', 6, 'LineWidth', 1.5, 'DisplayName', 'FEM Solution (Nodes)');
    
    % Add segment boundaries
    for i = 1:2
        x_pos = i * L;
        line([x_pos, x_pos], ylim, 'Color', [0.5, 0.5, 0.5], 'LineStyle', '--', 'LineWidth', 1);
        text(x_pos + 15, max(displacement_analytical) * 0.9, sprintf('Segment %d/%d boundary', i, i+1), ...
            'Rotation', 90, 'VerticalAlignment', 'top', 'FontSize', 9);
    end
    
    % Annotate with key values
    text(L/2, displacement_analytical(round(num_points/6)) * 1.1, ...
        sprintf('Max displacement in Segment 1: %.2f mm', max(displacement_analytical(x_analytical <= L))), ...
        'FontSize', 9, 'HorizontalAlignment', 'center', ...
        'BackgroundColor', [0.9, 0.9, 1], 'EdgeColor', 'blue', 'Margin', 3);
    
    text(L + L/2, displacement_analytical(round(num_points/2)) * 1.1, ...
        sprintf('Max displacement in Segment 2: %.2f mm', max(displacement_analytical(x_analytical <= 2*L & x_analytical > L))), ...
        'FontSize', 9, 'HorizontalAlignment', 'center', ...
        'BackgroundColor', [0.9, 0.9, 1], 'EdgeColor', 'blue', 'Margin', 3);
    
    text(2*L + L/2, displacement_analytical(round(5*num_points/6)) * 1.1, ...
        sprintf('Max displacement in Segment 3: %.2f mm', max(displacement_analytical(x_analytical > 2*L))), ...
        'FontSize', 9, 'HorizontalAlignment', 'center', ...
        'BackgroundColor', [0.9, 0.9, 1], 'EdgeColor', 'blue', 'Margin', 3);
    
    xlabel('Position along bar (mm)', 'FontSize', 12);
    ylabel('Displacement (mm)', 'FontSize', 12);
    title('Displacement Field Comparison: Analytical vs. FEM', 'FontSize', 14);
    legend('Location', 'northwest');
    grid on;
    
    % Save the displacement plot
    saveas(gcf, fullfile(plots_dir, 'displacement_field.png'));
    
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
    
    % Add stress annotations
    text(L/2, stress_analytical(round(num_points/6)) * 0.7, ...
        sprintf('Stress range in Segment 1: %.0f - %.0f MPa', ...
        min(stress_analytical(x_analytical <= L)), max(stress_analytical(x_analytical <= L))), ...
        'FontSize', 9, 'HorizontalAlignment', 'center', ...
        'BackgroundColor', [1, 0.9, 0.9], 'EdgeColor', 'red', 'Margin', 3);
    
    text(L + L/2, stress_analytical(round(num_points/2)) * 0.7, ...
        sprintf('Stress in Segment 2: %.0f MPa', stress_analytical(round(num_points/2))), ...
        'FontSize', 9, 'HorizontalAlignment', 'center', ...
        'BackgroundColor', [1, 0.9, 0.9], 'EdgeColor', 'red', 'Margin', 3);
    
    text(2*L + L/2, stress_analytical(round(5*num_points/6)) * 0.7, ...
        sprintf('Stress in Segment 3: %.0f MPa', stress_analytical(round(5*num_points/6))), ...
        'FontSize', 9, 'HorizontalAlignment', 'center', ...
        'BackgroundColor', [1, 0.9, 0.9], 'EdgeColor', 'red', 'Margin', 3);
    
    xlabel('Position along bar (mm)', 'FontSize', 12);
    ylabel('Stress (MPa)', 'FontSize', 12);
    title('Stress Field Comparison: Analytical vs. FEM', 'FontSize', 14);
    legend('Location', 'northeast');
    grid on;
    
    % Save the stress plot
    saveas(gcf, fullfile(plots_dir, 'stress_field.png'));
    
    %% CONVERGENCE STUDY
    % Perform convergence study with different mesh sizes
    element_counts = [4, 8, 16, 32];
    errors = zeros(size(element_counts));
    
    figure('Position', [100, 100, 900, 600]);
    
    for i = 1:length(element_counts)
        n_elements = element_counts(i);
        
        % Solve with current mesh size
        [elements_conv, nodes_conv, displacements_conv, stresses_conv] = ...
            solve_fem(A1, A2, A3, E1, E2, L, F1, F2, F3, n_elements);
        
        % Calculate element center positions
        x_centers_conv = zeros(length(elements_conv), 1);
        for j = 1:length(elements_conv)
            node1 = elements_conv(j, 1);
            node2 = elements_conv(j, 2);
            x_centers_conv(j) = (nodes_conv(node1) + nodes_conv(node2)) / 2;
        end
        
        % Calculate exact stresses at element centers
        exact_stresses = zeros(size(x_centers_conv));
        for j = 1:length(x_centers_conv)
            exact_stresses(j) = interp1(x_analytical, stress_analytical, x_centers_conv(j));
        end
        
        % Calculate RMS error
        errors(i) = sqrt(mean((stresses_conv - exact_stresses).^2));
        
        % Plot stress for this mesh size
        subplot(2, 2, i);
        plot(x_analytical, stress_analytical, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Analytical');
        hold on;
        plot(x_centers_conv, stresses_conv, 'ro', 'MarkerSize', 4, 'DisplayName', 'FEM');
        title(sprintf('%d Elements per Segment (RMS Error: %.2f MPa)', n_elements, errors(i)), 'FontSize', 10);
        xlabel('Position (mm)', 'FontSize', 9);
        ylabel('Stress (MPa)', 'FontSize', 9);
        legend('Location', 'northeast', 'FontSize', 8);
        grid on;
    end
    
    % Save the convergence plot
    saveas(gcf, fullfile(plots_dir, 'convergence_study.png'));
    
    % Plot error vs. element count
    figure('Position', [100, 100, 600, 400]);
    loglog(element_counts, errors, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
    
    % Add linear fit in log-log space to determine convergence rate
    p = polyfit(log(element_counts), log(errors), 1);
    convergence_rate = p(1);
    fitted_errors = exp(polyval(p, log(element_counts)));
    
    hold on;
    loglog(element_counts, fitted_errors, 'r--', 'LineWidth', 1.5);
    text(element_counts(2), fitted_errors(2)*1.2, sprintf('Convergence Rate: %.2f', convergence_rate), ...
        'FontSize', 10, 'FontWeight', 'bold');
    
    xlabel('Number of Elements per Segment', 'FontSize', 12);
    ylabel('RMS Error in Stress (MPa)', 'FontSize', 12);
    title('Convergence Study: Error vs. Mesh Refinement', 'FontSize', 14);
    grid on;
    
    % Save the error convergence plot
    saveas(gcf, fullfile(plots_dir, 'error_convergence.png'));
    
    %% CROSS-SECTIONAL AREA VISUALIZATION
    x_plot = linspace(0, total_length, 100);
    area_values = zeros(size(x_plot));
    
    % Calculate area at each point
    for i = 1:length(x_plot)
        area_values(i) = get_area_at_x(x_plot(i), A1, A2, A3, L);
    end
    
    figure('Position', [100, 100, 900, 400]);
    plot(x_plot, area_values, 'b-', 'LineWidth', 3);
    
    % Add segment boundaries
    for i = 1:2
        x_pos = i * L;
        line([x_pos, x_pos], ylim, 'Color', [0.5, 0.5, 0.5], 'LineStyle', '--', 'LineWidth', 1);
        text(x_pos + 15, mean(ylim), sprintf('Segment %d/%d boundary', i, i+1), ...
            'Rotation', 90, 'VerticalAlignment', 'middle', 'FontSize', 9);
    end
    
    % Annotate area values
    text(L/2, A1*0.85, sprintf('A1 to A2: %.0f to %.0f mm²', A1, A2), ...
        'FontSize', 10, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    text(L + L/2, A2*0.85, sprintf('A2: %.0f mm²', A2), ...
        'FontSize', 10, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    text(2*L + L/2, A3*0.85, sprintf('A3: %.0f mm²', A3), ...
        'FontSize', 10, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    
    xlabel('Position along bar (mm)', 'FontSize', 12);
    ylabel('Cross-sectional Area (mm²)', 'FontSize', 12);
    title('Cross-sectional Area Distribution', 'FontSize', 14);
    grid on;
    
    % Save the area plot
    saveas(gcf, fullfile(plots_dir, 'cross_sectional_area.png'));
    
    %% COMBINED VISUALIZATION
    figure('Position', [100, 100, 1000, 800]);
    
    % Subplot 1: Cross-sectional area
    subplot(3, 1, 1);
    plot(x_plot, area_values, 'b-', 'LineWidth', 2);
    for i = 1:2
        x_pos = i * L;
        line([x_pos, x_pos], ylim, 'Color', [0.5, 0.5, 0.5], 'LineStyle', '--', 'LineWidth', 1);
    end
    ylabel('Area (mm²)', 'FontSize', 10);
    title('Cross-sectional Area', 'FontSize', 12);
    grid on;
    
    % Subplot 2: Stress distribution
    subplot(3, 1, 2);
    plot(x_analytical, stress_analytical, 'r-', 'LineWidth', 2);
    hold on;
    plot(x_element_centers, element_stresses, 'ko', 'MarkerSize', 4);
    for i = 1:2
        x_pos = i * L;
        line([x_pos, x_pos], ylim, 'Color', [0.5, 0.5, 0.5], 'LineStyle', '--', 'LineWidth', 1);
    end
    ylabel('Stress (MPa)', 'FontSize', 10);
    title('Stress Distribution', 'FontSize', 12);
    legend('Analytical', 'FEM', 'Location', 'northeast', 'FontSize', 8);
    grid on;
    
    % Subplot 3: Displacement field
    subplot(3, 1, 3);
    plot(x_analytical, displacement_analytical, 'g-', 'LineWidth', 2);
    hold on;
    plot(x_fem, u_fem, 'ko', 'MarkerSize', 4);
    for i = 1:2
        x_pos = i * L;
        line([x_pos, x_pos], ylim, 'Color', [0.5, 0.5, 0.5], 'LineStyle', '--', 'LineWidth', 1);
    end
    xlabel('Position along bar (mm)', 'FontSize', 10);
    ylabel('Displacement (mm)', 'FontSize', 10);
    title('Displacement Field', 'FontSize', 12);
    grid on;
    
    % Overall title
    sgtitle('Composite Bar Structure: Complete Analysis', 'FontSize', 14, 'FontWeight', 'bold');
    
    % Save the combined visualization
    saveas(gcf, fullfile(plots_dir, 'combined_visualization.png'));
end

% Function to calculate analytical solution
function [displacement, stress] = analytical_solution(A1, A2, A3, E1, E2, L, F1, F2, F3, x)
    % Calculate analytical solution for displacement and stress
    displacement = zeros(size(x));
    stress = zeros(size(x));
    
    % Total length
    total_length = 3 * L;
    
    % Convert to compatible units
    E1 = E1 * 1000;  % GPa to MPa
    E2 = E2 * 1000;  % GPa to MPa
    F1 = F1 * 1000;  % kN to N
    F2 = F2 * 1000;  % kN to N
    F3 = F3 * 1000;  % kN to N
    
    % Calculate u at x = L (junction of segments 1 and 2)
    P1 = F1 + F2 + F3;  % Force in segment 1
    
    % Calculate u at x = 2L (junction of segments 2 and 3)
    P2 = F2 + F3;  % Force in segment 2
    
    % Calculate u at x = 3L (end of bar)
    P3 = F3;  % Force in segment 3
    
    % Integration constants
    % Find C1 from boundary condition u(0) = 0
    C1 = 0;
    
    % For segment 1: 0 <= x <= L
    % Linearly varying area: A(x) = A1 - (A1-A2)*x/L
    % Calculate u1(L) for use in determining C2
    u1_at_L = P1 * L / (E1 * A2) * log(A1/A2);
    
    % For segment 2: L <= x <= 2L
    % Constant area A2
    % Calculate u2(L) and use continuity u1(L) = u2(L) to find C2
    C2 = u1_at_L;
    
    % Calculate u2(2L) for use in determining C3
    u2_at_2L = C2 + P2 * L / (E1 * A2);
    
    % For segment 3: 2L <= x <= 3L
    % Constant area A3 and different modulus E2
    % Calculate u3(2L) and use continuity u2(2L) = u3(2L) to find C3
    C3 = u2_at_2L;
    
    % Calculate displacement and stress at each position
    for i = 1:length(x)
        xi = x(i);
        
        if xi <= L
            % Segment 1: linearly varying area
            A_x = A1 - (A1 - A2) * (xi / L);
            
            % Displacement
            displacement(i) = P1 * L / (E1 * A2) * log(A1 / A_x);
            
            % Stress
            stress(i) = P1 / A_x;
            
        elseif xi <= 2*L
            % Segment 2: constant area A2
            displacement(i) = C2 + P2 * (xi - L) / (E1 * A2);
            stress(i) = P2 / A2;
            
        else
            % Segment 3: constant area A3, different modulus E2
            displacement(i) = C3 + P3 * (xi - 2*L) / (E2 * A3);
            stress(i) = P3 / A3;
        end
    end
end

% Function to get area at position x
function area = get_area_at_x(x, A1, A2, A3, L)
    total_length = 3 * L;
    
    if x < 0 || x > total_length
        error('x must be between 0 and %f', total_length);
    end
    
    if x <= L
        % First segment: linearly varying area from A1 to A2
        area = A1 - (A1 - A2) * (x / L);
    elseif x <= 2 * L
        % Second segment: constant area A2
        area = A2;
    else
        % Third segment: constant area A3
        area = A3;
    end
end

% Function to solve using FEM
function [elements, nodes, nodal_displacements, element_stresses] = solve_fem(A1, A2, A3, E1, E2, L, F1, F2, F3, num_elements_per_segment)
    % Solve using Finite Element Method
    
    % Convert to compatible units
    E1 = E1 * 1000;  % GPa to MPa
    E2 = E2 * 1000;  % GPa to MPa
    F1 = F1 * 1000;  % kN to N
    F2 = F2 * 1000;  % kN to N
    F3 = F3 * 1000;  % kN to N
    
    % Total length and number of elements
    total_length = 3 * L;
    num_total_elements = 3 * num_elements_per_segment;
    num_nodes = num_total_elements + 1;
    
    % Node coordinates
    nodes = linspace(0, total_length, num_nodes);
    
    % Element connectivity (node1, node2)
    elements = zeros(num_total_elements, 2);
    for i = 1:num_total_elements
        elements(i, :) = [i, i+1];
    end
    
    % Initialize global stiffness matrix and force vector
    K = zeros(num_nodes, num_nodes);
    F = zeros(num_nodes, 1);
    
    % Assemble global stiffness matrix
    for e = 1:num_total_elements
        % Get element nodes
        node1 = elements(e, 1);
        node2 = elements(e, 2);
        
        % Element length
        x1 = nodes(node1);
        x2 = nodes(node2);
        Le = x2 - x1;
        
        % Element center for material property calculation
        x_center = (x1 + x2) / 2;
        
        % Determine element properties based on segment
        if x_center <= L
            % Segment 1: linearly varying area
            Ae = get_area_at_x(x_center, A1, A2, A3, L);
            Ee = E1;
        elseif x_center <= 2*L
            % Segment 2: constant area A2
            Ae = A2;
            Ee = E1;
        else
            % Segment 3: constant area A3, different modulus E2
            Ae = A3;
            Ee = E2;
        end
        
        % Element stiffness matrix
        ke = Ae * Ee / Le * [1, -1; -1, 1];
        
        % Assemble into global stiffness matrix
        K(node1, node1) = K(node1, node1) + ke(1, 1);
        K(node1, node2) = K(node1, node2) + ke(1, 2);
        K(node2, node1) = K(node2, node1) + ke(2, 1);
        K(node2, node2) = K(node2, node2) + ke(2, 2);
    end
    
    % Apply external forces at appropriate nodes
    % Find nodes closest to force application points
    [~, node_F1] = min(abs(nodes - L));
    [~, node_F2] = min(abs(nodes - 2*L));
    [~, node_F3] = min(abs(nodes - 3*L));
    
    F(node_F1) = F(node_F1) + F1;
    F(node_F2) = F(node_F2) + F2;
    F(node_F3) = F(node_F3) + F3;
    
    % Apply boundary condition: fixed at x = 0
    % Modify K and F for the boundary condition
    K(1, :) = 0;
    K(1, 1) = 1;
    F(1) = 0;
    
    % Solve for displacements
    nodal_displacements = K \ F;
    
    % Calculate element stresses
    element_stresses = zeros(num_total_elements, 1);
    for e = 1:num_total_elements
        % Get element nodes
        node1 = elements(e, 1);
        node2 = elements(e, 2);
        
        % Element length
        x1 = nodes(node1);
        x2 = nodes(node2);
        Le = x2 - x1;
        
        % Element center for material property calculation
        x_center = (x1 + x2) / 2;
        
        % Determine element properties based on segment
        if x_center <= L
            % Segment 1: linearly varying area
            Ae = get_area_at_x(x_center, A1, A2, A3, L);
            Ee = E1;
        elseif x_center <= 2*L
            % Segment 2: constant area A2
            Ae = A2;
            Ee = E1;
        else
            % Segment 3: constant area A3, different modulus E2
            Ae = A3;
            Ee = E2;
        end
        
        % Element strain: B*u
        B = [-1/Le, 1/Le];
        u_element = [nodal_displacements(node1); nodal_displacements(node2)];
        strain = B * u_element;
        
        % Element stress: E*strain
        element_stresses(e) = Ee * strain;
    end
end
