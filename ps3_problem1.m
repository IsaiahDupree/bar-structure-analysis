%% FEA Solver for Problem Set 3 - Problem 1: Three 2-node Elements
% This script solves Problem 1 from Problem Set 3

%% Clear workspace and command window
clear all;
close all;
clc;

%% Create output directory for figures
if ~exist('figures', 'dir')
    mkdir('figures');
end

%% Problem parameters
node_positions = [0, 2, 4, 6];  % x-coordinates of nodes
displacements = [1, 3, 4, 6];   % Given displacements at nodes

fprintf('Solving Problem 1: Three 2-node Elements\n');

%% 1. Derive shape functions for all three elements
fprintf('\n1. Shape Functions for 2-node elements:\n');
fprintf('N₁(ξ) = (1-ξ)/2\n');
fprintf('N₂(ξ) = (1+ξ)/2\n');

% For reference, ξ goes from -1 to 1 in each element
% Element 1: nodes 1 and 2 (x = 0 to 2)
% Element 2: nodes 2 and 3 (x = 2 to 4)
% Element 3: nodes 3 and 4 (x = 4 to 6)

% Element lengths
L1 = node_positions(2) - node_positions(1);  % 2
L2 = node_positions(3) - node_positions(2);  % 2
L3 = node_positions(4) - node_positions(3);  % 2

%% 2. Find displacements in each element at ξ = 0.5
xi = 0.5;  % Given ξ value

% Shape function values at ξ = 0.5
N1_at_xi = (1 - xi) / 2;  % 0.25
N2_at_xi = (1 + xi) / 2;  % 0.75

% Displacements at ξ = 0.5 in each element
u_elem1_at_xi = N1_at_xi * displacements(1) + N2_at_xi * displacements(2);
u_elem2_at_xi = N1_at_xi * displacements(2) + N2_at_xi * displacements(3);
u_elem3_at_xi = N1_at_xi * displacements(3) + N2_at_xi * displacements(4);

fprintf('\n2. Displacements at ξ = 0.5:\n');
fprintf('Element 1: u(ξ=0.5) = %.4f\n', u_elem1_at_xi);
fprintf('Element 2: u(ξ=0.5) = %.4f\n', u_elem2_at_xi);
fprintf('Element 3: u(ξ=0.5) = %.4f\n', u_elem3_at_xi);

%% 3. Find strains in each element at ξ = 0.5
% For a 2-node element, strain is constant: ε = (u₂-u₁)/L
% dN₁/dξ = -0.5, dN₂/dξ = 0.5 (derivatives of shape functions)
% Jacobian for each element: J = L/2
% dN₁/dx = dN₁/dξ * dξ/dx = -0.5 * 2/L = -1/L
% dN₂/dx = dN₂/dξ * dξ/dx = 0.5 * 2/L = 1/L

strain_elem1 = (displacements(2) - displacements(1)) / L1;
strain_elem2 = (displacements(3) - displacements(2)) / L2;
strain_elem3 = (displacements(4) - displacements(3)) / L3;

fprintf('\n3. Strains at ξ = 0.5:\n');
fprintf('Element 1: ε(ξ=0.5) = %.4f\n', strain_elem1);
fprintf('Element 2: ε(ξ=0.5) = %.4f\n', strain_elem2);
fprintf('Element 3: ε(ξ=0.5) = %.4f\n', strain_elem3);

%% 4. Calculate the stiffness matrix of all elements using Gauss integration
% For simplicity, we'll assume E (Young's modulus) = 1 and A (Cross-sectional area) = 1
E = 1.0;  % Young's modulus (assumed)
A = 1.0;  % Cross-sectional area (assumed)

% Calculate stiffness matrices
k1 = calculate_stiffness_matrix(L1, E, A);
k2 = calculate_stiffness_matrix(L2, E, A);
k3 = calculate_stiffness_matrix(L3, E, A);

fprintf('\n4. Element Stiffness Matrices:\n');
fprintf('Element 1 Stiffness Matrix:\n');
disp(k1);
fprintf('Element 2 Stiffness Matrix:\n');
disp(k2);
fprintf('Element 3 Stiffness Matrix:\n');
disp(k3);

% Generate plots
generate_plots_problem1(node_positions, displacements, u_elem1_at_xi, u_elem2_at_xi, ...
                      u_elem3_at_xi, strain_elem1, strain_elem2, strain_elem3, k1, k2, k3);

%% Helper function to calculate stiffness matrix
function k = calculate_stiffness_matrix(L, E, A)
    % For 2-node elements, the stiffness matrix is:
    % k = [E*A/L, -E*A/L; -E*A/L, E*A/L]
    k = E * A / L * [1, -1; -1, 1];
end

%% Helper function to generate plots
function generate_plots_problem1(node_positions, displacements, u_elem1_at_xi, u_elem2_at_xi, ...
                               u_elem3_at_xi, strain_elem1, strain_elem2, strain_elem3, k1, k2, k3)
    % Plot 1: Bar configuration and displacements
    figure('Position', [100, 100, 800, 600]);
    
    % x-coordinates for plotting
    x_plot = linspace(0, node_positions(end), 300);
    
    % Displacement function
    u_plot = zeros(size(x_plot));
    for i = 1:length(x_plot)
        x = x_plot(i);
        if x <= node_positions(2)  % Element 1
            xi = 2 * (x - node_positions(1)) / (node_positions(2) - node_positions(1)) - 1;
            N1 = (1 - xi) / 2;
            N2 = (1 + xi) / 2;
            u_plot(i) = N1 * displacements(1) + N2 * displacements(2);
        elseif x <= node_positions(3)  % Element 2
            xi = 2 * (x - node_positions(2)) / (node_positions(3) - node_positions(2)) - 1;
            N1 = (1 - xi) / 2;
            N2 = (1 + xi) / 2;
            u_plot(i) = N1 * displacements(2) + N2 * displacements(3);
        else  % Element 3
            xi = 2 * (x - node_positions(3)) / (node_positions(4) - node_positions(3)) - 1;
            N1 = (1 - xi) / 2;
            N2 = (1 + xi) / 2;
            u_plot(i) = N1 * displacements(3) + N2 * displacements(4);
        end
    end
    
    % Plot displacement field
    plot(x_plot, u_plot, 'b-', 'LineWidth', 2);
    hold on;
    plot(node_positions, displacements, 'ro', 'MarkerSize', 8);
    
    % Mark ξ = 0.5 points on each element
    x_elem1_mid = node_positions(1) + 0.75 * (node_positions(2) - node_positions(1));
    x_elem2_mid = node_positions(2) + 0.75 * (node_positions(3) - node_positions(2));
    x_elem3_mid = node_positions(3) + 0.75 * (node_positions(4) - node_positions(3));
    
    plot(x_elem1_mid, u_elem1_at_xi, 'g*', 'MarkerSize', 10);
    plot(x_elem2_mid, u_elem2_at_xi, 'g*', 'MarkerSize', 10);
    plot(x_elem3_mid, u_elem3_at_xi, 'g*', 'MarkerSize', 10);
    
    text(x_elem1_mid, u_elem1_at_xi+0.2, ['\xi=0.5, u=', num2str(u_elem1_at_xi, '%.4f')], ...
        'HorizontalAlignment', 'center', 'FontSize', 9);
    text(x_elem2_mid, u_elem2_at_xi+0.2, ['\xi=0.5, u=', num2str(u_elem2_at_xi, '%.4f')], ...
        'HorizontalAlignment', 'center', 'FontSize', 9);
    text(x_elem3_mid, u_elem3_at_xi+0.2, ['\xi=0.5, u=', num2str(u_elem3_at_xi, '%.4f')], ...
        'HorizontalAlignment', 'center', 'FontSize', 9);
    
    % Add node annotations
    for i = 1:length(node_positions)
        text(node_positions(i), displacements(i)-0.5, ['Node ', num2str(i), sprintf('\nu=%g', displacements(i))], ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
            'BackgroundColor', 'white', 'EdgeColor', 'gray', 'Margin', 3);
    end
    
    % Add element divisions
    for i = 2:length(node_positions)-1
        line([node_positions(i) node_positions(i)], ylim, 'Color', 'gray', 'LineStyle', '--', 'LineWidth', 1);
    end
    
    % Add element labels
    for i = 1:length(node_positions)-1
        mid_x = (node_positions(i) + node_positions(i+1)) / 2;
        text(mid_x, max(displacements)+0.7, ['Element ', num2str(i)], 'HorizontalAlignment', 'center', ...
            'BackgroundColor', 'lightyellow', 'EdgeColor', 'gray', 'Margin', 3);
    end
    
    xlabel('Position x');
    ylabel('Displacement u');
    title('Problem 1: Displacement Field for Three 2-node Elements');
    grid on;
    grid minor;
    legend('Displacement Field', 'Nodal Displacements', 'Location', 'best');
    
    % Save figure
    saveas(gcf, fullfile('figures', 'ps3_problem1_displacement.png'));
    
    % Plot 2: Element strains
    figure('Position', [100, 100, 600, 400]);
    
    elements = {'Element 1', 'Element 2', 'Element 3'};
    strains = [strain_elem1, strain_elem2, strain_elem3];
    
    bar(1:3, strains);
    colormap('summer');
    
    % Add strain values
    for i = 1:length(strains)
        text(i, strains(i)+0.02, num2str(strains(i), '%.4f'), 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom');
    end
    
    set(gca, 'XTick', 1:3, 'XTickLabel', elements);
    ylabel('Strain');
    title('Problem 1: Element Strains at \xi = 0.5');
    grid on;
    grid minor;
    
    % Save figure
    saveas(gcf, fullfile('figures', 'ps3_problem1_strains.png'));
end
