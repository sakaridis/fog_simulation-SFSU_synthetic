function t = disparity_to_transmittance_homogeneous_cityscapes(...
    disparity_in_pixels_complete, camera_parameters_filename, beta)
%DISPARITY_TO_TRANSMITTANCE_HOMOGENEOUS_CITYSCAPES  Convert disparity map of
%Cityscapes image (in pixels) that is free from holes to transmittance map for
%the given attenuation coefficient.
%
%   INPUTS:
%
%   -|disparity_in_pixels_complete|: matrix with elements of type double, with
%    the same resolution as Cityscapes RGB images. Its elements are positive
%    numbers.
%
%   -|camera_parameters_file|: full path to JSON file where camera parameters
%    are stored.
%
%   -|beta|: attenuation coefficient. Positive scalar.
%
%   OUTPUTS:
%
%   -|t|: matrix of same size and type as |disparity_in_pixels_complete|,
%    containing the values of transmittance at each pixel, which range inside
%    [0, 1].

% Add the required paths.
current_function_full_name = mfilename('fullpath');
addpath(fullfile(fileparts(current_function_full_name), '..', 'utilities'));
addpath_relative_to_caller(current_function_full_name,...
    fullfile('..', 'Fog_simulation'));

t = transmission_exponential(distance_in_meters_cityscapes(...
    depth_in_meters_cityscapes_disparity_in_pixels(...
    disparity_in_pixels_complete, camera_parameters_filename),...
    camera_parameters_filename), beta);

end

