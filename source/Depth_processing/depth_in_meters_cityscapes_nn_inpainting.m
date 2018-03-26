function depth_map_in_meters =...
    depth_in_meters_cityscapes_nn_inpainting(input_disparity,...
    camera_parameters_file, varargin)
%DEPTH_IN_METERS_CITYSCAPES_NN_INPAINTING  Compute complete depth map in meters
%from provided Cityscapes disparity map and relevant camera parameters, using
%nearest neighbor (NN) inpainting of the missing parts from the known ones.
%
%   INPUTS:
%
%   -|input_disparity|: matrix of uint16 format with same resolution as
%   Cityscapes RGB images.
%
%   -|camera_parameters_file|: full path to JSON file where camera parameters
%   are stored.
%
%   OUTPUTS:
%
%   -|depth_map_in_meters|: matrix of double format containing depth in meters.

% Input raw depth map provided in Cityscapes.
[depth_input, is_depth_input_invalid] =...
    depth_in_meters_cityscapes_with_invalid_parts(input_disparity,...
    camera_parameters_file);

% Find the pixels with known depth that SURROUND the invalid parts: they are the
% only candidates for nearest neighbors of the invalid pixels out of the entire
% set of pixels with known depth.
is_depth_input_invalid_dilated =...
    imdilate(is_depth_input_invalid, strel('square', 3));
is_depth_input_invalid_outer_boundary =...
    is_depth_input_invalid_dilated & ~is_depth_input_invalid;
inds_boundary = find(is_depth_input_invalid_outer_boundary);

% Get the coordinates of "boundary" and invalid pixels for nearest neighbor
% search.
[height, width] = size(input_disparity);
[X, Y] = meshgrid(1:width, 1:height);
pixel_coords = [X(:), Y(:)];
pixel_coords_boundary = pixel_coords(inds_boundary, :);
pixel_coords_invalid = pixel_coords(is_depth_input_invalid, :);

% Find nearest neighbors for all invalid pixels.
inds_nearest_known_pixels =...
    knnsearch(pixel_coords_boundary, pixel_coords_invalid);

% Inpaint invalid parts using the values at the nearest neighbor pixels.
depth_map_in_meters = depth_input;
depth_map_in_meters(is_depth_input_invalid) =...
    depth_input(inds_boundary(inds_nearest_known_pixels));

end

