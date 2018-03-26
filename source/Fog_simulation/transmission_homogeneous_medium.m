function t = transmission_homogeneous_medium(d, beta, camera_parameters_file)
%TRANSMISSION_HOMOGENEOUS_MEDIUM  Compute transmission map using given depth
%map, based on the Beer-Lambert law. Distinguish between scene depth |d| and
%distance between camera and object depicted at each pixel, |l|.
%   INPUTS:
%
%   -|d|: H-by-W matrix with values of depth for processed image in meters.
%
%   -|beta|: attenuation coefficient in inverse meters. Constant, since the
%    medium is homogeneous.
%
%   OUTPUTS:
%
%   -|t|: H-by-W matrix with medium transmission values ranging in [0, 1].

% Add directory with |distance_in_meters_cityscapes| function to path.
current_script_full_name = mfilename('fullpath');
current_script_directory = fileparts(current_script_full_name);
addpath(fullfile(current_script_directory, '..', 'Depth_processing'));

% Compute scene distance from camera for each pixel.
l = distance_in_meters_cityscapes(d, camera_parameters_file);

% Beer-Lambert law for homogeneous medium.
t = exp(-beta * l);

end

