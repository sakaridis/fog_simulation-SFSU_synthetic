function [depth_map_in_meters, is_depth_invalid] =...
    depth_in_meters_cityscapes_stereoscopic_inpainting(input_disparity,...
    camera_parameters_file, left_image, left_image_uint8, right_image, varargin)
%DEPTH_IN_METERS_CITYSCAPES_STEREOSCOPIC_INPAINTING  Compute denoised and
%complete depth map in meters from provided Cityscapes disparity map, left and
%right RGB image of stereo pair and relevant camera parameters, using the method
%proposed in the Stereoscopic Inpainting paper combined with SLIC for
%segmentation and a novel objective for segment matching.
%
%   INPUTS:
%
%   -|input_disparity|: matrix of uint16 format with same resolution as
%   Cityscapes RGB images.
%
%   -|camera_parameters_file|: full path to JSON file where camera parameters
%   are stored.
%
%   -|left_image|: left image of stereo pair, the depth of which is being
%   processed, in double format.
%
%   -|left_image_uint8|: left image of stereo pair, the depth of which is being
%   processed, in uint8 format.
%
%   -|right_image|: right image of stereo pair in double format. This image is
%   auxiliary.
%
%   OUTPUTS:
%
%   -|depth_map_in_meters|: matrix of double format containing depth in meters.
%   It may contain erroneous parts that are known to be such: they correspond to
%   SLIC segments with bad result from Stereoscopic Inpainting and assume zero
%   values at some parts of the segments.
%
%   -|is_depth_invalid|: boolean matrix, where true indicates that the
%   corresponding pixel in |depth_map_in_meters| certainly contains a erroneous
%   value.

% Suppress the nearly singular matrix warning, which otherwise occurs during
% RANSAC plane fitting.
warning_state = warning;
warning('off', 'MATLAB:nearlySingularMatrix');

% Input the required data.
[depth_input, is_depth_input_invalid] =...
    depth_in_meters_cityscapes_with_invalid_parts(input_disparity,...
    camera_parameters_file);
left_disparity = disparity_in_pixels_cityscapes(input_disparity);
[height, width, ~] = size(left_image);
[X, Y] = meshgrid(1:width, 1:height);
points = [X(:), Y(:)];


% Photo-consistency check to create binary mask |M_L| with unreliable depth.
epsilon = 12 / 255;
photoconsistency_outlier_mask = outliers_photoconsistency(left_image,...
    right_image, left_disparity, epsilon);
M_L = photoconsistency_outlier_mask | is_depth_input_invalid;


% Segmentation of left image in superpixel-wise segments using SLIC.
% 1) Set parameters.
number_of_preferred_segments = 2048;
m = 10;
% 2) Perform segmentation. |S_L| is the image containing the segmentation result
%    for the left image. |K| is the number of output segments.
[S_L, K] = slicmex(left_image_uint8, number_of_preferred_segments, m);
% 3) Convert segmentation result to double and one-based format for subsequent
%    operations.
S_L = double(S_L) + 1;


% Segment classification and plane fitting with RANSAC for adequately visible
% segments.
minimum_count_known = 20;
minimum_fraction_known = 0.6;
visible = false(1, K);
% Threshold for large-margin outliers w.r.t. depth.
thresh_large = 50;
% RANSAC parameters.
alpha = 10 ^ -2;
iter = 2000;
p = 1 - 10 ^ -2;
% Fix random seed for reproducibility.
rng('default');
% Initialize output depth map to the input depth read from Cityscapes.
depth_map_in_meters = depth_input;
% Initialize matrix with plane parameters.
planes = zeros(3, K);
% Initialize boolean vector denoting whether the plane of segment lies
% infinitely deep.
is_plane_at_inf = false(1, K);
% Convert input image to CIELAB colorspace for subsequent color similarity
% computation.
left_image_lab = rgb2lab(left_image);
left_image_l = left_image_lab(:, :, 1);
left_image_a = left_image_lab(:, :, 2);
left_image_b = left_image_lab(:, :, 3);
% Initialize matrix with average CIELAB colors per segment.
lab_averages = zeros(K, 3);
% Initialize matrix with segment centroids.
centroids = zeros(K, 2);
% Main loop over all segments.
for i = 1:K
    
    % Select current segment.
    segment_current = S_L == i;
    
    % Identify the following disjoint subsets of current segment:
    % 1) Set of pixels with a known depth value.
    % 2) Set of pixels with an unreliable depth value.
    known_pixels = segment_current & ~M_L;
    count_pixels = nnz(segment_current);
    unreliable_pixels = segment_current & M_L;
    count_unreliable = nnz(unreliable_pixels);
    count_known = count_pixels - count_unreliable;
    
    if count_known >= max(minimum_count_known,...
            minimum_fraction_known * count_pixels)
        
        % Segment is adequately visible. Run RANSAC for plane fitting.
        
        visible(i) = 1;
        
        % Treat pixels with infinite depth separately. They are either all
        % inliers or all outliers and any inlier set can't contain both pixels
        % with infinite depth and others with finite depth. This observation
        % simplifies the application of RANSAC to the complete set of pixels
        % with a known depth value.
        known_pixels_finite = known_pixels & depth_input < Inf;
        count_finite = nnz(known_pixels_finite);
        count_inf = count_known - count_finite;
        
        % If the pixels with infinite depth are more than those with finite
        % depth, then they directly constitute the largest inlier set of the
        % segment.
        if count_inf > count_finite
            is_plane_at_inf(i) = true;
            planes(:, i) = Inf * ones(3, 1);
        else
            % Otherwise, apply RANSAC only to the pixels with finite depth and
            % check whether the resulting inlier set is larger than that from
            % infinite depth.
            depth_known_finite = depth_input(known_pixels_finite);
            
            [plane_finite, inlier_count_finite] =...
                ransac_depth_plane(points(known_pixels_finite, :),...
                depth_known_finite, alpha, iter, p);
            
            if count_inf > inlier_count_finite
                is_plane_at_inf(i) = true;
                planes(:, i) = Inf * ones(3, 1);
            else
                planes(:, i) = plane_finite;
            end
        end
        
        % Inpaint unreliable depth values using the fitted plane.
        depth_map_in_meters(unreliable_pixels) =...
            inpaint_depth_with_plane(points(unreliable_pixels, :),...
            planes(:, i));

        % Inpaint initially known depth values which strongly disagree with the
        % depth that is implied by the segment's plane, using the latter.
        depth_to_inpaint = inpaint_depth_with_plane(points(known_pixels, :),...
            planes(:, i));
        indices_known_pixels = find(known_pixels);
        large_outlier_pixels = abs(depth_to_inpaint -...
            depth_map_in_meters(known_pixels)) > thresh_large;
        depth_map_in_meters(indices_known_pixels(large_outlier_pixels)) =...
            depth_to_inpaint(large_outlier_pixels);
    end
    
    % Retrieve CIELAB values for all pixels of current segment.
    segment_lab = [left_image_l(segment_current),...
        left_image_a(segment_current), left_image_b(segment_current)];
    
    % Compute average CIELAB color, which is a meaningful operation as this
    % colorspace is perceptually uniform and therefore it can be treated as
    % Euclidean.
    lab_averages(i, :) = mean(segment_lab);
    
    % Retrieve pixel coordinates for all pixels of current segment and compute
    % the segment centroid.
    segment_xy = points(segment_current, :);
    centroids(i, :) = mean(segment_xy);
    
end
% Number of adequately visible segments.
visible_indices = find(visible);
V = length(visible_indices);
% Number of unreliable segments.
invisible_indices = find(~visible);
U = length(invisible_indices);


% Assignment of unreliable segments to visible ones with greedy matching.
% Color-similarity scale parameter.
lambda = m;
% Parameter to scale spatial distances between segment centroids. Approximately
% equal to the initial sampling step for superpixel centers of SLIC.
S = sqrt(height * width / double(K));
% Segregate invisible from visible segments following the Stereoscopic
% Inpainting paper.
lab_averages_invisible = lab_averages(invisible_indices, :);
lab_averages_visible = lab_averages(visible_indices, :);
centroids_invisible = centroids(invisible_indices, :);
centroids_visible = centroids(visible_indices, :);
% Pairwise squared distances in range domain.
lab_average_sqdists =...
    squared_euclidean_distances_exact([lab_averages_invisible;
    lab_averages_visible], lab_averages_invisible);
% Pairwise squared distances in spatial domain.
centroid_sqdists = squared_euclidean_distances_exact([centroids_invisible;
    centroids_visible], centroids_invisible);
% Form the complete objective as a weighted combination of the color similarity
% term and the spatial proximity term, where weighting is performed using the
% scale parameters for the two domains. Diagonal elements are set to infinity,
% in order to avoid matching invisible segments to themselves.
E = lab_average_sqdists + (lambda / S) ^ 2 * centroid_sqdists;
E(1:(U + V + 1):(U + V) * U) = Inf;
% Initialization of greedy algorithm for matching.
[E_min_vis, ind] = min(E(U + 1:end, :), [], 1);
E_min_invis = Inf(1, U);
E_min = E_min_vis;
unmatched = true(1, U);
best_match_is_visible = true(1, U);
ind_invisible = ind;
ind_final = zeros(size(ind));
is_matched_with_visible = true(1, U);
% Main loop for matching invisible segments to segments with already assigned
% planes.
while any(unmatched)
    
    % Identify unmatched segment with lowest value of objective.
    unmatched_indices = find(unmatched);
    [~, j_tmp] = min(E_min(unmatched));
    j = unmatched_indices(j_tmp);
    
    segment_current_id = invisible_indices(j);
    segment_current = S_L == segment_current_id;
    unreliable_pixels = segment_current & M_L;
    known_pixels = segment_current & ~M_L;
    
    % Assign an initially visible or invisible segment to the former segment,
    % based on the "source" of the objective.
    if best_match_is_visible(j)
        ind_final(j) = ind(j);
        planes(:, segment_current_id) =...
            planes(:, visible_indices(ind_final(j)));
        is_plane_at_inf(segment_current_id) =...
            is_plane_at_inf(visible_indices(ind_final(j)));
    else
        ind_final(j) = ind_invisible(j);
        planes(:, segment_current_id) =...
            planes(:, invisible_indices(ind_final(j)));
        is_plane_at_inf(segment_current_id) =...
            is_plane_at_inf(invisible_indices(ind_final(j)));
        is_matched_with_visible(j) = false;
    end
    
    % Inpaint unreliable depth values using the assigned plane.
    depth_map_in_meters(unreliable_pixels) =...
        inpaint_depth_with_plane(points(unreliable_pixels, :),...
        planes(:, segment_current_id));
    
    % Inpaint initially known depth values which strongly disagree with the
    % depth that is implied by the segment's plane, using the latter.
    depth_to_inpaint = inpaint_depth_with_plane(points(known_pixels, :),...
        planes(:, segment_current_id));
    indices_known_pixels = find(known_pixels);
    large_outlier_pixels =...
        abs(depth_to_inpaint - depth_map_in_meters(known_pixels)) >...
        thresh_large;
    depth_map_in_meters(indices_known_pixels(large_outlier_pixels)) =...
        depth_to_inpaint(large_outlier_pixels);
    
    % Update loop variables.
    ind_invisible(E_min_invis > E(j, :)) = j;
    E_min_invis = min(E_min_invis, E(j, :));
    E_min = min(E_min, E_min_invis);
    best_match_is_visible = E_min_vis <= E_min;
    unmatched(j) = false;
    
end


% Identify segments for which the filled depth values are erroneous.
% If a segment contains some not strictly positive depth values, depth
% completion has not worked correctly for that segment.
wrong_segments_ids = unique(S_L(depth_map_in_meters <= 0));
is_depth_invalid = ismember(S_L, wrong_segments_ids);

% Restore initial warning settings.
warning(warning_state);

end

