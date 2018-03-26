function [plane, inlier_count, inliers] = ransac_depth_plane(x, depth, alpha,...
    iter, p)
%RANSAC_DEPTH_PLANE  Fit a 3D depth plane to a set of input points using RANSAC.
%The distance of a point from the plane is defined as the vertical offset.
%
%   INPUTS:
%
%   -|x|: |P|-by-2 matrix with image-plane coordinates of input points.
%
%   -|depth|: |P|-by-1 vector with depth values for input points.
%
%   -|alpha|: positive constant regulating the inlier threshold used in RANSAC.
%
%   -|iter|: maximum number of RANSAC iterations.
%
%   -|p|: target probability for having obtained a pure inlier set, used to cut
%   down on number of RANSAC iterations.
%
%   OUTPUTS:
%
%   -|plane|: 3-by-1 vector holding the parameters of the fitted plane, so that
%   for each point [x, y, d] it holds d = [x, y, 1] * plane.
%
%   -|inlier_count|: cardinality of the final inlier set.
%
%   -|inliers|: |inlier_count|-by-1 vector with indices of the members of the
%   final inlier set. 

% Determine inlier threshold |theta| based on the depth values. Goal: invariance
% to range of depth values.
theta = alpha * median(depth);

% Number of plane parameters.
N = 3;

% Total number of points.
P = length(depth);

% Represent pixel coordinates homogeneously.
x_h = [x, ones(P, 1)];

% Initialize variables to save the best fit across iterations.
inliers = [];
inlier_count = 0;
inlier_ratio = 0;

% Main RANSAC loop.
for i = 1:iter
    colinear = true;
    while colinear
        % Select a sample at random.
        permut = randperm(P);
        sample = permut(1:N);
        x_h_sample = x_h(sample, :);
        
        % Proceed to plane fitting only if the points of the sample are not
        % colinear, i.e. square matrix with homogeneous coordinates of
        % points in the sample is of full rank.
        if rank(x_h_sample) == N
            colinear = false;
        end
    end
    
    % Fit a plane to the current sample.
    plane_tmp = fit_depth_plane(x_h_sample, depth(sample));
    
    % Identify inliers for current fit.
    
    % Vertical distance.
    inliers_tmp = abs(x_h * plane_tmp - depth) <= theta;
    % Orthogonal distance.
    % plane_tmp_normal = [plane_tmp(1:N-1); 1];
    % inliers_tmp = abs([points(known_pixels_finite, :), depth_known_finite] *...
    %     plane_tmp_normal - plane_tmp(N)) / norm(plane_tmp_normal) <= theta;
    inlier_tmp_count = nnz(inliers_tmp);
    
    % Update best fit.
    if inlier_tmp_count > inlier_count
        inlier_count = inlier_tmp_count;
        inlier_ratio = inlier_count / P;
        inliers = find(inliers_tmp);
    end
    
    % An adaptive version of RANSAC is implemented. If a pure inlier sample
    % has been already examined with probability greater that |p|, avoid
    % performing more iterations.
    p_bound = 1 - (1 - inlier_ratio ^ N) ^ i;
    if p <= p_bound
        break;
    end
end

% Perform a final least squares fit to compute the best depth plane with respect
% to all points in the maximum inlier set that has been found with RANSAC.
plane = fit_depth_plane(x_h(inliers, :), depth(inliers));

end

