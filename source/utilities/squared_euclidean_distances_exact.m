function squared_distances = squared_euclidean_distances_exact(points_1,...
    points_2)
%SQUARED_EUCLIDEAN_DISTANCES_EXACT  Compute squared Euclidean distance between
%each pair from two sets of points. Exact implementation.
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
points_1_rep = zeros(n, p, d);
points_2_rep = zeros(n, p, d);

for i = 1:d
    [points_2_rep(:, :, i), points_1_rep(:, :, i)] =...
        meshgrid(points_2(:, i), points_1(:, i));
end

squared_distances = sum((points_1_rep - points_2_rep) .^ 2, 3);

end

