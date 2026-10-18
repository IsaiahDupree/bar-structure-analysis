function h = arrow(p0,p1,varargin)
%ARROW Draw a line with an arrowhead.
%   ARROW(P0,P1) draws a line with an arrow from point P0 to point P1.
%   The P0 and P1 are points in the form [X,Y].
%
%   ARROW(P0,P1,PARAM,VALUE,...) draws the arrow with additional
%   parameters specified by PARAM,VALUE pairs:
%       'Length'    Length of the arrowhead (default: 10)
%       'Width'     Width of the arrowhead (default: 8)
%       'FaceColor' Color of the arrowhead (default: 'black')
%       'EdgeColor' Color of the arrow edges (default: 'black')
%       'LineWidth' Width of the arrow line (default: 1.5)
%
%   H = ARROW(...) returns a handle to the created objects.

% Define default parameters
params = struct(...
    'Length', 10, ...
    'Width', 8, ...
    'FaceColor', 'black', ...
    'EdgeColor', 'black', ...
    'LineWidth', 1.5);

% Process parameter-value pairs
for i = 1:2:length(varargin)
    if isfield(params, varargin{i})
        params.(varargin{i}) = varargin{i+1};
    end
end

% Draw the line
hLine = line([p0(1) p1(1)], [p0(2) p1(2)], 'Color', params.EdgeColor, 'LineWidth', params.LineWidth);

% Calculate the direction vector
v = p1 - p0;
v = v / norm(v);  % Normalize

% Calculate perpendicular vector
w = [-v(2), v(1)];

% Calculate arrowhead points
tip = p1;
baseLeft = p1 - params.Length * v + (params.Width/2) * w;
baseRight = p1 - params.Length * v - (params.Width/2) * w;

% Draw the arrowhead
hArrow = patch('XData', [tip(1), baseLeft(1), baseRight(1)], ...
               'YData', [tip(2), baseLeft(2), baseRight(2)], ...
               'FaceColor', params.FaceColor, ...
               'EdgeColor', params.EdgeColor, ...
               'LineWidth', params.LineWidth);

% Return handles if requested
if nargout > 0
    h = [hLine, hArrow];
end
end
