%% FEA Solver for Problem Set 3 - Problem 2: Two 3-node Elements
% This script solves Problem 2 from Problem Set 3

%% Clear workspace and command window
clear all;
close all;
clc;

%% Create output directory for figures
if ~exist('figures', 'dir')
    mkdir('figures');
end

%% Problem parameters
node_positions = [0, 2, 4, 6, 8];  % x-coordinates of nodes
displacements = [0, -1, 2, -1, 4];   % Given displacements at nodes

fprintf('Solving Problem 2: Two 3-node Elements\n');

%% 1. Derive shape functions for 3-node elements
fprintf('\n1. Shape Functions for 3-node elements:\n');
fprintf('N₁(ξ) = ξ(ξ-1)/2\n');
fprintf('N₂(ξ) = (1+ξ)(1-ξ)\n');
fprintf('N₃(ξ) = ξ(ξ+1)/2\n');

% For reference, ξ goes from -1 to 1 in each element
% Element 1: nodes 1, 2, 3 (x = 0, 2, 4) → ξ = -1, 0, 1
% Element 2: nodes 3, 4, 5 (x = 4, 6, 8) → ξ = -1, 0, 1

% Element lengths
L1 = node_positions(3) - node_positions(1);  % 4
L2 = node_positions(5) - node_positions(3);  % 4

%% 2. Find displacements in each element at ξ = -0.5
xi = -0.5;  % Given ξ value

% Shape function values at ξ = -0.5
N1_at_xi = xi * (xi - 1) / 2;  % 0.125
N2_at_xi = (1 + xi) * (1 - xi);  % 0.75
N3_at_xi = xi * (xi + 1) / 2;  % -0.125

% Displacements at ξ = -0.5 in each element
u_elem1_at_xi = N1_at_xi * displacements(1) + N2_at_xi * displacements(2) + N3_at_xi * displacements(3);
u_elem2_at_xi = N1_at_xi * displacements(3) + N2_at_xi * displacements(4) + N3_at_xi * displacements(5);

fprintf('\n2. Displacements at ξ = -0.5:\n');
fprintf('Element 1: u(ξ=-0.5) = %.4f\n', u_elem1_at_xi);
fprintf('Element 2: u(ξ=-0.5) = %.4f\n', u_elem2_at_xi);

%% 3. Find strains in each element at ξ = -0.5
% For a 3-node element:
% dN₁/dξ = ξ - 0.5
% dN₂/dξ = -2ξ
% dN₃/dξ = ξ + 0.5

% Derivatives of shape functions at ξ = -0.5
dN1_dxi = xi - 0.5;  % -1.0
dN2_dxi = -2 * xi;  % 1.0
dN3_dxi = xi + 0.5;  % 0.0

% Jacobian for element 1: J = L1/2 = 2
% Jacobian for element 2: J = L2/2 = 2
J1 = L1 / 2;
J2 = L2 / 2;

% Strain = B * u where B = [dN₁/dx dN₂/dx dN₃/dx]
% dN/dx = dN/dξ * dξ/dx = dN/dξ * 1/J

% Element 1 strain at ξ = -0.5
B1 = [dN1_dxi, dN2_dxi, dN3_dxi] / J1;
strain_elem1 = B1 * displacements(1:3)';

% Element 2 strain at ξ = -0.5
B2 = [dN1_dxi, dN2_dxi, dN3_dxi] / J2;
strain_elem2 = B2 * displacements(3:5)';

fprintf('\n3. Strains at ξ = -0.5:\n');
fprintf('Element 1: ε(ξ=-0.5) = %.4f\n', strain_elem1);
fprintf('Element 2: ε(ξ=-0.5) = %.4f\n', strain_elem2);

%% 4. Calculate the stiffness matrix of all elements using integration points ξ = -0.5 and ξ = 0.5
% For simplicity, we'll assume E (Young's modulus) = 1 and A (Cross-sectional area) = 1
E = 1.0;  % Young's modulus (assumed)
A = 1.0;  % Cross-sectional area (assumed)

% Integration points
xi_points = [-0.5, 0.5];
weights = [1.0, 1.0];

% Calculate stiffness matrices
k1 = calculate_stiffness_matrix(L1, E, A, xi_points, weights);
k2 = calculate_stiffness_matrix(L2, E, A, xi_points, weights);

fprintf('\n4. Element Stiffness Matrices:\n');
fprintf('Element 1 Stiffness Matrix:\n');
disp(k1);
fprintf('Element 2 Stiffness Matrix:\n');
disp(k2);

% Generate plots
generate_plots_problem2(node_positions, displacements, u_elem1_at_xi, u_elem2_at_xi, ...
                      strain_elem1, strain_elem2, k1, k2);

%% Helper function to calculate stiffness matrix
function k = calculate_stiffness_matrix(L, E, A, xi_points, weights)
    % Calculate stiffness matrix for a 3-node element using numerical integration
    k = zeros(3, 3);
    
    % Jacobian: J = L/2
    J = L / 2;
    
    for i = 1:length(xi_points)
        xi = xi_points(i);
        % Derivatives of shape functions at this integration point
        dN1_dxi = xi - 0.5;
        dN2_dxi = -2 * xi;
        dN3_dxi = xi + 0.5;
        
        % B matrix at this point: B = [dN₁/dx dN₂/dx dN₃/dx]
        B = [dN1_dxi, dN2_dxi, dN3_dxi] / J;
        
        % Contribution to stiffness matrix: w * B^T * E * B * A * det(J)
        k = k + weights(i) * B' * E * B * A * J;
    end
end

%% Helper function to generate plots
function generate_plots_problem2(node_positions, displacements, u_elem1_at_xi, u_elem2_at_xi, ...
                               strain_elem1, strain_elem2, k1, k2)
    % Plot 1: Bar configuration and displacements
    figure('Position', [100, 100, 900, 600]);
    
    % x-coordinates for plotting
    x_plot = linspace(0, node_positions(end), 300);
    
    % Displacement function
    u_plot = zeros(size(x_plot));
    for i = 1:length(x_plot)
        x = x_plot(i);
        if x <= node_positions(3)  % Element 1
            % Map x to ξ in the element
            xi = 2 * (x - node_positions(1)) / (node_positions(3) - node_positions(1)) - 1;
            
            % Compute shape functions
            N1 = xi * (xi - 1) / 2;
            N2 = (1 + xi) * (1 - xi);
            N3 = xi * (xi + 1) / 2;
            
            u_plot(i) = N1 * displacements(1) + N2 * displacements(2) + N3 * displacements(3);
            
        else  % Element 2
            % Map x to ξ in the element
            xi = 2 * (x - node_positions(3)) / (node_positions(5) - node_positions(3)) - 1;
            
            % Compute shape functions
            N1 = xi * (xi - 1) / 2;
            N2 = (1 + xi) * (1 - xi);
            N3 = xi * (xi + 1) / 2;
            
            u_plot(i) = N1 * displacements(3) + N2 * displacements(4) + N3 * displacements(5);
        end
    end
    
    % Plot displacement field
    plot(x_plot, u_plot, 'b-', 'LineWidth', 2);
    hold on;
    plot(node_positions, displacements, 'ro', 'MarkerSize', 8);
    
    % Mark ξ = -0.5 points on each element
    x_elem1_mid = node_positions(1) + 0.25 * (node_positions(3) - node_positions(1));  % x at ξ = -0.5 in element 1
    x_elem2_mid = node_positions(3) + 0.25 * (node_positions(5) - node_positions(3));  % x at ξ = -0.5 in element 2
    
    plot(x_elem1_mid, u_elem1_at_xi, 'g*', 'MarkerSize', 10);
    plot(x_elem2_mid, u_elem2_at_xi, 'g*', 'MarkerSize', 10);
    
    text(x_elem1_mid, u_elem1_at_xi+0.3, ['\xi=-0.5, u=', num2str(u_elem1_at_xi, '%.4f')], ...
        'HorizontalAlignment', 'center', 'FontSize', 9);
    text(x_elem2_mid, u_elem2_at_xi+0.3, ['\xi=-0.5, u=', num2str(u_elem2_at_xi, '%.4f')], ...
        'HorizontalAlignment', 'center', 'FontSize', 9);
    
    % Add node annotations
    for i = 1:length(node_positions)
        text(node_positions(i), displacements(i)+0.5*(-1)^(i+1), ...
            ['Node ', num2str(i), sprintf('\nu=%g', displacements(i))], ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'center', ...
            'BackgroundColor', 'white', 'EdgeColor', 'gray', 'Margin', 3);
    end
    
    % Add element divisions
    line([node_positions(3) node_positions(3)], ylim, 'Color', 'gray', 'LineStyle', '--', 'LineWidth', 1);
    
    % Add element labels
    mid_x1 = (node_positions(1) + node_positions(3)) / 2;
    mid_x2 = (node_positions(3) + node_positions(5)) / 2;
    
    text(mid_x1, max(displacements)+1.0, 'Element 1 (3-node)', 'HorizontalAlignment', 'center', ...
        'BackgroundColor', 'lightyellow', 'EdgeColor', 'gray', 'Margin', 3);
    text(mid_x2, max(displacements)+1.0, 'Element 2 (3-node)', 'HorizontalAlignment', 'center', ...
        'BackgroundColor', 'lightyellow', 'EdgeColor', 'gray', 'Margin', 3);
    
    xlabel('Position x');
    ylabel('Displacement u');
    title('Problem 2: Displacement Field for Two 3-node Elements');
    grid on;
    grid minor;
    legend('Displacement Field', 'Nodal Displacements', 'Location', 'best');
    
    % Set y-axis limits with some padding
    ylim([min(displacements)-1.0, max(displacements)+1.5]);
    
    % Save figure
    saveas(gcf, fullfile('figures', 'ps3_problem2_displacement.png'));
    
    % Plot 2: Element strains
    figure('Position', [100, 100, 600, 400]);
    
    elements = {'Element 1', 'Element 2'};
    strains = [strain_elem1, strain_elem2];
    
    bar(1:2, strains);
    colormap('summer');
    
    % Add strain values
    for i = 1:length(strains)
        text(i, strains(i)+0.05*sign(strains(i)), num2str(strains(i), '%.4f'), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom');
    end
    
    set(gca, 'XTick', 1:2, 'XTickLabel', elements);
    ylabel('Strain');
    title('Problem 2: Element Strains at \xi = -0.5');
    grid on;
    grid minor;
    
    % Save figure
    saveas(gcf, fullfile('figures', 'ps3_problem2_strains.png'));
    
    % Plot 3: Element stiffness matrices visualization as heatmaps
    figure('Position', [100, 100, 900, 400]);
    
    % Plot stiffness matrix for element 1
    subplot(1, 2, 1);
    imagesc(k1);
    title('Element 1 Stiffness Matrix');
    colorbar;
    axis equal;
    axis tight;
    set(gca, 'XTick', 1:3, 'YTick', 1:3);
    set(gca, 'XTickLabel', {'Node 1', 'Node 2', 'Node 3'});
    set(gca, 'YTickLabel', {'Node 1', 'Node 2', 'Node 3'});
    colormap('jet');
    
    % Add text annotations to show matrix values
    for i = 1:3
        for j = 1:3
            text(j, i, num2str(k1(i,j), '%.3f'), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                'Color', 'white');
        end
    end
    
    % Plot stiffness matrix for element 2
    subplot(1, 2, 2);
    imagesc(k2);
    title('Element 2 Stiffness Matrix');
    colorbar;
    axis equal;
    axis tight;
    set(gca, 'XTick', 1:3, 'YTick', 1:3);
    set(gca, 'XTickLabel', {'Node 3', 'Node 4', 'Node 5'});
    set(gca, 'YTickLabel', {'Node 3', 'Node 4', 'Node 5'});
    colormap('jet');
    
    % Add text annotations to show matrix values
    for i = 1:3
        for j = 1:3
            text(j, i, num2str(k2(i,j), '%.3f'), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                'Color', 'white');
        end
    end
    
    % Save figure
    saveas(gcf, fullfile('figures', 'ps3_problem2_stiffness.png'));
end
