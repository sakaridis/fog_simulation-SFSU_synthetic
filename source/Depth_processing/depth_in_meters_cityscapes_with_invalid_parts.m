function [depth_map_in_meters, is_disparity_invalid] =...
    depth_in_meters_cityscapes_with_invalid_parts(input_disparity,...
    camera_parameters_file)
%DEPTH_IN_METERS_CITYSCAPES_WITH_INVALID_PARTS  Compute depth map in meters from
%provided Cityscapes disparity map and relevant camera parameters.
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
%   It contains invalid measurements that are known to be such and are set to
%   infinity by convention, and may also contain more erroneous values that are
%   not known to be such and can assume any positive value.
%
%   -|is_disparity_invalid|: boolean matrix, where true indicates that the
%   corresponding pixel in |depth_map_in_meters| contains a wrong infinity
%   value.

% Identify known wrong values for disparity, based on specifications provided in
% the Cityscapes README.
is_disparity_invalid = input_disparity == 0;

% Compute the disparity in pixels, based on specifications provided in the
% Cityscapes README.
disparity_in_pixels = disparity_in_pixels_cityscapes(input_disparity);
is_disparity_zero = disparity_in_pixels == 0;

% Retrieve baseline and focal length in x-axis, which are both required in order
% to get the absolute scale of the depth map.
[B, f_x] = camera_parameters_cityscapes(camera_parameters_file);

% Compute the depth as inversely proportional to disparity. Wherever the
% disparity is zero, the depth is equal to infinity. Convention: any known
% invalid disparity is assigned infinite depth.
depth_map_in_meters = zeros(size(disparity_in_pixels));
depth_map_in_meters(~is_disparity_invalid & ~is_disparity_zero) =...
    B * f_x ./ disparity_in_pixels(~is_disparity_invalid & ~is_disparity_zero);
depth_map_in_meters(is_disparity_zero | is_disparity_invalid) = Inf;

end

