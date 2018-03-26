function squared_distances = squared_euclidean_distances_fast(points_1,...
    points_2)
%SQUARED_EUCLIDEAN_DISTANCES_FAST  Compute squared Euclidean distance between
%each pair from two sets of points. Fast implementation.
%
%   INPUTS:
%
%   -|points_1|: |n|-by-|d| matrix with first set of points.
%
%   -|points_2|: |p|-by-|d| matrix with second set of points.
%
%   OUTPUTS:
%
%   -|squared_distances|: |n|-by-|p| matrix with squared Euclidean distances
%   between each pair of points from the two sets.

[n, d] = size(points_1);
p = size(points_2, 1);

% Compute squared distance for every pair of points.
tmp_1 = zeros(n, 3 * d);
tmp_2 = zeros(p, 3 * d);
for i = 1:d
    tmp_1(:, 3 * i - 2:3 * i) = [ones(n, 1), -2 * points_1(:, i),...
        points_1(:, i) .^ 2];
    tmp_2(:, 3 * i - 2:3 * i) = [points_2(:, i) .^ 2 ,...
        points_2(:, i), ones(p, 1)];
end
squared_distances = tmp_1 * tmp_2.';

end

