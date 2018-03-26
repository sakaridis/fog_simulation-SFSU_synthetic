function depth_map_in_meters =...
    depth_in_meters_cityscapes_disparity_in_pixels(disparity_in_pixels,...
    camera_parameters_file)
%DEPTH_IN_METERS_CITYSCAPES_DISPARITY_IN_PIXELS  Compute depth map in meters
%from complete (preprocessed) disparity map in pixels and relevant camera
%parameters.
%
%   INPUTS:
%
%   -|disparity_in_pixels|: matrix with elements of type double, with the same
%    resolution as Cityscapes RGB images. Its elements are positive numbers.
%
%   -|camera_parameters_file|: full path to JSON file where camera parameters
%    are stored.
%
%   OUTPUTS:
%
%   -|depth_map_in_meters|: matrix of double format containing depth in meters.
%    Its elements are equal to |Inf| for corresponding elements in
%    |disparity_in_pixels| that are equal to 0.

is_disparity_zero = disparity_in_pixels == 0;

% Retrieve baseline and focal length in x-axis, which are both required in order
% to get the absolute scale of the depth map.
[B, f_x] = camera_parameters_cityscapes(camera_parameters_file);

% Compute the depth as inversely proportional to disparity.
depth_map_in_meters = zeros(size(disparity_in_pixels));
depth_map_in_meters(~is_disparity_zero) =...
    B * f_x ./ disparity_in_pixels(~is_disparity_zero);
depth_map_in_meters(is_disparity_zero) = Inf;

end

