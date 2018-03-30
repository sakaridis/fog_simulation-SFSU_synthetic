%% Input the required data

close all;
clear;

% Add required paths.

current_script_full_name = mfilename('fullpath');
current_script_directory = fileparts(current_script_full_name);
addpath(fullfile(current_script_directory, '..', 'utilities'));
addpath_relative_to_caller(current_script_full_name,...
    fullfile('..', 'Fog_simulation'));
addpath_relative_to_caller(current_script_full_name,...
    fullfile('..', 'Depth_processing'));
addpath_relative_to_caller(current_script_full_name,...
    fullfile('..', 'Dehazing', 'Dark_channel_prior'));
addpath_relative_to_caller(current_script_full_name,...
    fullfile('..', 'external', 'SLIC_mex'));

% Flags for plotting figures.
plot_basic_figures = 1;
plot_many_figures = 0;

mask_color = [1, 0, 0];

images_root_dir = fullfile(current_script_directory, '..', '..', 'data',...
    'demos');

% Example images. Uncomment whichever you would like to experiment with.
image_basename = 'hamburg_000000_046078';
% image_basename = 'dusseldorf_000020_000019';
% image_basename = 'cologne_000088_000019';
% image_basename = 'aachen_000001_000019';

left_input = imread(fullfile(images_root_dir, 'leftImg8bit',...
    strcat(image_basename, '_leftImg8bit.png')));
right_input = imread(fullfile(images_root_dir, 'rightImg8bit',...
    strcat(image_basename, '_rightImg8bit.png')));
left_input_disparity = imread(fullfile(images_root_dir, 'disparity',...
    strcat(image_basename, '_disparity.png')));
camera_parameters_file = fullfile(images_root_dir, 'camera',...
    strcat(image_basename, '_camera.json'));

left_image = im2double(left_input);
right_image = im2double(right_input);
[height, width, ~] = size(left_image);
[X, Y] = meshgrid(1:width, 1:height);
points = [X(:), Y(:)];

[depth, is_depth_invalid] =...
    depth_in_meters_cityscapes_with_invalid_parts(left_input_disparity,...
    camera_parameters_file);
left_disparity = disparity_in_pixels_cityscapes(left_input_disparity);

%% Photo-consistency check to create binary mask |M_L| with unreliable depth

epsilon = 12 / 255;
photoconsistency_outlier_mask = outliers_photoconsistency(left_image,...
    right_image, left_disparity, epsilon);

M_L = photoconsistency_outlier_mask | is_depth_invalid;

% Optionally plot mask of unreliable depth pixels.
if plot_basic_figures
    figure;
    imshow(M_L);
end

%% Segmentation of left image in superpixel-wise segments

% SLIC segmentation:

% 1) Set parameters.
number_of_preferred_segments = 2048;
m = 10;

% 2) Perform segmentation. |S_L| is the image containing the segmentation result
%    for the left image. |K| is the number of output segments.
[S_L, K] = slicmex(left_input, number_of_preferred_segments, m);

% Optionally plot the result.
if plot_basic_figures
    S_L_help = imfilter(S_L, fspecial('average'), 'replicate');
    superpixel_boundaries = S_L ~= S_L_help;
    left_image_superpixel_boundaries = mask_color_image(left_image,...
        superpixel_boundaries, mask_color);
    figure;
    imshow(left_image_superpixel_boundaries);
end

% Convert segmentation result to double and one-based format for subsequent
% operations.
S_L = double(S_L) + 1;

% Convert input image to L*a*b* colorspace for subsequent color similarity
% computation.
MATLAB_release = version('-release');
first_release_with_rgb2lab = '2014b';
if cstrcmp(MATLAB_release, first_release_with_rgb2lab) >= 0
    left_image_lab = rgb2lab(left_image);
else
    cform = makecform('srgb2lab');
    left_image_lab = applycform(left_image, cform);
end
left_image_l = left_image_lab(:, :, 1);
left_image_a = left_image_lab(:, :, 2);
left_image_b = left_image_lab(:, :, 3);

%% Segment classification and plane fitting with RANSAC for reliable segments

minimum_count_known = 20;
minimum_fraction_known = 0.6;
visible = false(1, K);
completely_visible = false(1, K);

% Threshold for large-margin outliers w.r.t. depth.
thresh_large = 50;

% RANSAC parameters.
alpha = 10 ^ -2;
iter = 2000;
p = 1 - 10 ^ -2;

% Suppress the nearly singular matrix warning, which otherwise occurs during
% RANSAC plane fitting.
warning_state = warning;
warning('off', 'MATLAB:nearlySingularMatrix');

% Fix random seed for reproducibility.
rng('default');

% Initialize output depth map to the input depth read from Cityscapes.
depth_filled = depth;

% Monitor total number of depth values that are filled in the plane fitting
% step.
pixels_filled_fitting = 0;

% Initialize matrix with plane parameters.
planes = zeros(3, K);
% Initialize boolean vector denoting whether the plane of segment lies
% infinitely deep.
is_plane_at_inf = false(1, K);

% Initialize matrix with average L*a*b* colors per segment.
lab_averages = zeros(K, 3);

% Initialize matrix with segment centroids.
centroids = zeros(K, 2);

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
    
    if count_unreliable == 0
        completely_visible(i) = 1;
    end
    
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
        known_pixels_finite = known_pixels & depth < Inf;
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
            depth_known_finite = depth(known_pixels_finite);
            
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
        depth_filled(unreliable_pixels) =...
            inpaint_depth_with_plane(points(unreliable_pixels, :),...
            planes(:, i));
        
        pixels_filled_fitting = pixels_filled_fitting + count_unreliable;
        
        if ~is_plane_at_inf(i) && ~mod(i, 200) && plot_many_figures
            figure;
            histogram(depth_known_finite);
            figure;
            histogram(depth_filled(unreliable_pixels));
        end

        % Inpaint initially known depth values which strongly disagree with the
        % depth that is implied by the segment's plane, using the latter.
        depth_to_inpaint = inpaint_depth_with_plane(points(known_pixels, :),...
            planes(:, i));
        indices_known_pixels = find(known_pixels);
        large_outlier_pixels = abs(depth_to_inpaint -...
            depth_filled(known_pixels)) > thresh_large;
        depth_filled(indices_known_pixels(large_outlier_pixels)) =...
            depth_to_inpaint(large_outlier_pixels);
    end
    
    % Retrieve L*a*b* values for all pixels of current segment.
    segment_lab = [left_image_l(segment_current),...
        left_image_a(segment_current), left_image_b(segment_current)];
    
    % Compute average L*a*b* color, which is a meaningful operation as this
    % colorspace is perceptually uniform and therefore it can be treated as
    % Euclidean.
    lab_averages(i, :) = mean(segment_lab);
    
    % Retrieve pixel coordinates for all pixels of current segment and compute
    % the segment centroid.
    segment_xy = points(segment_current, :);
    centroids(i, :) = mean(segment_xy);
    
end

% Restore initial warning settings.
warning(warning_state);

visible_indices = find(visible);
invisible_indices = find(~visible);
completely_visible_indices = find(completely_visible);

% Number of adequately visible segments.
V = length(visible_indices);

% Number of unreliable segments.
U = length(invisible_indices);

% Optionally plot union of invisible segments as a binary image.
if plot_basic_figures
    figure;
    imshow(ismember(S_L, invisible_indices));
end

%% Assignment of unreliable segments to reliable ones with greedy matching

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
matched_count = 0;

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
    depth_filled(unreliable_pixels) =...
        inpaint_depth_with_plane(points(unreliable_pixels, :),...
        planes(:, segment_current_id));
    
    % Inpaint initially known depth values which strongly disagree with the
    % depth that is implied by the segment's plane, using the latter.
    depth_to_inpaint = inpaint_depth_with_plane(points(known_pixels, :),...
        planes(:, segment_current_id));
    indices_known_pixels = find(known_pixels);
    large_outlier_pixels =...
        abs(depth_to_inpaint - depth_filled(known_pixels)) > thresh_large;
    depth_filled(indices_known_pixels(large_outlier_pixels)) =...
        depth_to_inpaint(large_outlier_pixels);
    
    % Update loop variables.
    ind_invisible(E_min_invis > E(j, :)) = j;
    E_min_invis = min(E_min_invis, E(j, :));
    E_min = min(E_min, E_min_invis);
    best_match_is_visible = E_min_vis <= E_min;
    unmatched(j) = false;
    
    matched_count = matched_count + 1;
end

%% Inspect how the depth values look like in a hole-free segment

segment_current_id = completely_visible_indices(20);
segment_current = S_L == segment_current_id;
segment_sample_translated_dr = imtranslate(segment_current, [1, 1]);
segment_sample_translated_r = imtranslate(segment_current, [1, 0]);
segment_sample_translated_d = imtranslate(segment_current, [0, 1]);
upper_right_triangles_origins = segment_current &...
    segment_sample_translated_dr & segment_sample_translated_r;
lower_left_triangles_origins = segment_current &...
    segment_sample_translated_dr & segment_sample_translated_d;
tri = [find(upper_right_triangles_origins),...
    find(imtranslate(upper_right_triangles_origins, [-1, -1]) &...
    segment_current),...
    find(imtranslate(upper_right_triangles_origins, [-1, 0]) & segment_current);
    find(lower_left_triangles_origins),...
    find(imtranslate(lower_left_triangles_origins, [-1, -1]) &...
    segment_current),...
    find(imtranslate(lower_left_triangles_origins, [0, -1]) & segment_current)];
standard_deviation_of_depth_inside_segment = std(depth(segment_current));

if plot_many_figures
    figure;
    trisurf(tri, points(:, 1), points(:, 2), depth(:), 'FaceColor', 'interp');
    figure;
    histogram(depth(segment_current));
end

%% Identify segments for which the filled depth values are erroneous

% If a segment contains some not strictly positive depth values, depth
% completion has not worked correctly for that segment.
wrong_segments_ids = unique(S_L(depth_filled <= 0));
is_depth_filled_invalid = ismember(S_L, wrong_segments_ids);

%% Fog simulation

% Transmittance map. Compute 3 different versions:
% 1) from raw Cityscapes depth. Holes are treated as infinitely deep parts.
% 2) from denoised and filled Cityscapes depth.
% 3) by applying guided filtering to 2).
beta = 0.01;

% 1)
t = transmission_homogeneous_medium(depth, beta, camera_parameters_file);

% 2)
t_filled = transmission_homogeneous_medium(depth_filled, beta,...
    camera_parameters_file);

% 3) Apply guided image filtering to the transmittance map that has been
%    computed from the filled depth map, using the original clean image as the
%    guidance image.
window_size = 41;
mu = 1e-3;
t_guided = clip_to_unit_range(transmission_guided_filtering(t_filled,...
    left_image, window_size, mu));

% Foggy image using transmittance map and atmospheric light. Try 3 different
% ways to set the value of atmospheric light:
% 1) manual setting.
% 2) estimation with method proposed in Dark Channel Prior paper.
% 3) estimation with modification of Dark Channel Prior method proposed in
%    Regression Forests paper.

% 1)
c = 0.96;
L = repmat(c, 1, 1, 3);

% 2)
neighborhood_size_dark_channel = 15;
left_image_dark = get_dark_channel(left_image, neighborhood_size_dark_channel);
L_dcp = estimate_atmospheric_light_dcp(left_image_dark, left_image);

% 3)
L_rf = estimate_atmospheric_light_rf(left_image_dark, left_image);

I = haze_linear(left_image, t, L);
I_filled = haze_linear(left_image, t_filled, L);
I_guided = haze_linear(left_image, t_guided, L);
I_guided_dcp = haze_linear(left_image, t_guided, L_dcp);
I_guided_rf = haze_linear(left_image, t_guided, L_rf);

%% Show fog synthesis result with marked wrong disparities

I_annotated_unreliable = mask_color_image(I, M_L, mask_color);
I_annotated_invalid = mask_color_image(I, is_depth_invalid, mask_color);
I_annotated_inconsistent_only = mask_color_image(I,...
    ~photoconsistency_outlier_mask & ~is_depth_invalid, mask_color);

if plot_basic_figures
    figure;
    imshow(t_guided);
    figure;
    imshow(I_guided_rf);
end

if plot_many_figures
    figure;
    imshow(t);
    figure;
    imshow(t_filled);
    figure;
    imshow(left_image);
    figure;
    imshow(I);
    figure;
    imshow(I_filled);
    figure;
    imshow(I_guided);
    figure;
    imshow(I_guided_dcp);
    figure;
    imshow(is_depth_filled_invalid);
    figure;
    imshowpair(I, I_filled, 'montage');
    figure;
    imshowpair(I, I_annotated_unreliable, 'montage');
    figure;
    imshowpair(I, I_annotated_invalid, 'montage');
    figure;
    imshowpair(I, I_annotated_inconsistent_only, 'montage');
end

