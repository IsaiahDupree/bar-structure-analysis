function postPlots(mesh, U, stress, P)
% POSTPLOTS Generate plots for displacement and stress
%
% Parameters:
% mesh - Matrix with columns [node1, node2, x0, length]
% U - Vector of nodal displacements
% stress - Vector of element stresses
% P - Structure containing problem parameters

% Create plots directory if it doesn't exist
if ~exist('plots', 'dir')
    mkdir('plots');
end

% Create x coordinates for plotting
x = linspace(0, 3*P.L, length(U));

% Figure 1: Nodal Displacement
figure('Position', [100, 100, 800, 500]);
plot(x, U, '-o', 'LineWidth', 1.5);
xlabel('x (mm)');
ylabel('u (mm)');
title('Nodal Displacement');
grid on;
saveas(gcf, 'plots/displacement_field.png');

% Figure 2: Element Stress
figure('Position', [100, 100, 800, 500]);
% Create x coordinates for stress visualization
xStress = [mesh(:,1); mesh(end,2)] * P.L / mesh(end,2);
stressPlot = [stress; stress(end)];  % Repeat the last stress value for proper stairs plot

stairs(xStress, stressPlot, 'LineWidth', 1.5);
xlabel('x (mm)');
ylabel('\sigma (MPa)');
title('Element Stress');
grid on;

% Add segment boundaries
hold on;
xline(P.L, 'r--', 'LineWidth', 1.5);
xline(2*P.L, 'r--', 'LineWidth', 1.5);
hold off;

saveas(gcf, 'plots/stress_field.png');

% Figure 3: Cross-sectional Area Distribution
figure('Position', [100, 100, 800, 500]);
x_plot = linspace(0, 3*P.L, 1000);
area = zeros(size(x_plot));

for i = 1:length(x_plot)
    xi = x_plot(i);
    if xi <= P.L  % First segment
        % Linear interpolation from A1 to A2
        area(i) = P.A(1) + (P.A(2) - P.A(1)) * (xi / P.L);
    elseif xi <= 2 * P.L  % Second segment
        area(i) = P.A(2);
    else  % Third segment
        area(i) = P.A(3);
    end
end

plot(x_plot, area, 'LineWidth', 1.5);
hold on;
xline(P.L, 'r--', 'LineWidth', 1.5);
xline(2*P.L, 'r--', 'LineWidth', 1.5);
title('Cross-sectional Area Distribution');
xlabel('x (mm)');
ylabel('Area (mm²)');
grid on;
saveas(gcf, 'plots/area_distribution.png');

fprintf('Plots generated successfully and saved to the "plots" directory.\n');

end
