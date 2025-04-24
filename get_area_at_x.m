function area = get_area_at_x(x, A1, A2, A3, L)
%GET_AREA_AT_X Get the cross-sectional area at position x
%   
%   For segment 1, the area varies linearly from A1 to A2
%   For segment 2, the area is constant A2
%   For segment 3, the area is constant A3
%
%   Parameters:
%   -----------
%   x : float
%       Position along the bar
%   A1, A2, A3 : float
%       Cross-sectional areas for each segment
%   L : float
%       Length of each segment
%
%   Returns:
%   --------
%   area : float
%       Cross-sectional area at position x

% Total length of the bar
total_length = 3 * L;

% Check if x is within valid range
if x < 0 || x > total_length
    error('x must be between 0 and %f', total_length);
end

% Calculate area based on position
if x <= L  % First segment
    % Linear interpolation from A1 to A2
    area = A1 + (A2 - A1) * (x / L);
elseif x <= 2 * L  % Second segment
    area = A2;
else  % Third segment
    area = A3;
end

end
