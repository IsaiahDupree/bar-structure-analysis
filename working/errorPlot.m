function errorPlot(mesh, sigmaFE, sigmaExactFn, P)
% ERRORPLOT Generate error plot comparing FEM and analytical solutions
%
% Parameters:
% mesh - Matrix with columns [node1, node2, x0, length]
% sigmaFE - Vector of element stresses from FEM solution
% sigmaExactFn - Function handle for analytical stress calculation
% P - Structure containing problem parameters

% Calculate element center coordinates
xCenter = ((mesh(:,1) + mesh(:,2)) / 2) * P.L / mesh(end,2);

% Calculate analytical stresses at element centers
sigmaAnalytical = sigmaExactFn(xCenter, P);

% Calculate relative error (percentage)
err = abs((sigmaFE - sigmaAnalytical) ./ sigmaAnalytical) * 100;

% Create error plot
figure('Position', [100, 100, 800, 500]);
stem(xCenter, err, 'filled');
ylabel('|\epsilon| (%)');
xlabel('x (mm)');
title('FE vs Analytical Stress Error');
grid on;

% Add horizontal line at 5% error threshold
hold on;
yline(5, 'r--', '5% Threshold', 'LineWidth', 1.5);
hold off;

% Save the plot
saveas(gcf, 'plots/error_plot.png');

% Print summary statistics
fprintf('Error Analysis:\n');
fprintf('  Maximum error: %.2f%%\n', max(err));
fprintf('  Average error: %.2f%%\n', mean(err));
fprintf('  Number of elements exceeding 5%% error: %d (out of %d)\n', ...
        sum(err > 5), length(err));

% Check if mesh refinement is needed
if max(err) > 5
    fprintf('  RECOMMENDATION: Refine mesh where error exceeds 5%%\n');
    
    % Identify segments needing refinement
    need_refinement = false(1,3);
    
    for s = 1:3
        % Calculate the element range for each segment
        if s == 1
            elem_range = 1:P.nEl(1);
        elseif s == 2
            elem_range = (P.nEl(1)+1):(P.nEl(1)+P.nEl(2));
        else
            elem_range = (P.nEl(1)+P.nEl(2)+1):(P.nEl(1)+P.nEl(2)+P.nEl(3));
        end
        
        % Check if any element in this segment exceeds 5% error
        if any(err(elem_range) > 5)
            need_refinement(s) = true;
            fprintf('    Segment %d needs refinement\n', s);
        end
    end
else
    fprintf('  Mesh accuracy is sufficient (all errors below 5%%)\n');
end

end

% Helper function: Analytical stress calculation
function s = analyticStress(x, P)
    % Calculate total reaction at fixed end
    R = sum(P.F);
    
    % Internal forces in each segment
    N1 = R;                    % Segment 1: N1 = R
    N2 = R - P.F(1);           % Segment 2: N2 = R - F1
    N3 = R - P.F(1) - P.F(2);  % Segment 3: N3 = R - F1 - F2

    % Calculate stress at each point
    s = zeros(size(x));
    for i = 1:numel(x)
        if x(i) <= P.L
            % Segment 1: Variable cross-section
            A = P.A(1) + (P.A(2) - P.A(1)) * (x(i) / P.L);
            s(i) = N1 / A;
        elseif x(i) <= 2*P.L
            % Segment 2: Constant cross-section A2
            s(i) = N2 / P.A(2);
        else
            % Segment 3: Constant cross-section A3
            s(i) = N3 / P.A(3);
        end
    end
end
