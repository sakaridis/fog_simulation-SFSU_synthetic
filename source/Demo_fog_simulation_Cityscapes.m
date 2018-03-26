% Demo of fog simulation pipeline presented in our article
% "Semantic Foggy Scene Understanding with Synthetic Data"
% (International Journal of Computer Vision, 2018)
% for application to real clear-weather outdoor scenes from the Cityscapes
% dataset.

clear;

% ------------------------------------------------------------------------------

% Add required paths.

current_script_full_name = mfilename('fullpath');
current_script_directory = fileparts(current_script_full_name);
addpath(fullfile(current_script_directory, 'utilities'));
addpath_relative_to_caller(current_script_full_name,...
    fullfile('Fog_simulation'));
addpath_relative_to_caller(current_script_full_name,...
    fullfile('Depth_processing'));
addpath_relative_to_caller(current_script_full_name,...
    fullfile('Dehazing', 'Dark_channel_prior'));
addpath_relative_to_caller(current_script_full_name,...
    fullfile('external', 'SLIC_mex'));

% ------------------------------------------------------------------------------

% Input the required data.

images_root_dir = fullfile(current_script_directory, '..', 'data', 'demos');

% Example images. Uncomment whichever image you like to experiment with.
image_basename = 'hamburg_000000_046078';
% image_basename = 'dusseldorf_000020_000019';
% image_basename = 'cologne_000088_000019';
% image_basename = 'aachen_000001_000019';

left_image_uint8 = imread(fullfile(images_root_dir, 'leftImg8bit',...
    strcat(image_basename, '_leftImg8bit.png')));
right_image = im2double(imread(fullfile(images_root_dir, 'rightImg8bit',...
    strcat(image_basename, '_rightImg8bit.png'))));
left_disparity = imread(fullfile(images_root_dir, 'disparity',...
    strcat(image_basename, '_disparity.png')));
camera_parameters_file = fullfile(images_root_dir, 'camera',...
    strcat(image_basename, '_camera.json'));

% Bring original, clear left image to double precision for subsequent
% computations.
left_image = im2double(left_image_uint8);

% ------------------------------------------------------------------------------

% Configure result files and directories, and parameters for fog simulation.

results_root_dir = fullfile(current_script_directory, '..', 'output', 'demos');
depth_results_dir = 'depth_stereoscopic';
foggy_results_dir = 'leftImg8bit_foggy';
transmittance_results_dir = 'leftImg8bit_transmittance';

% Attenuation coefficient (in inverse meters).
% In case of fog, the value of the attenuation coefficient is greater or equal
% to 0.003.
beta = 0.01;

suffix_depth = '_depth_stereoscopic';
suffix_foggy = strcat('_foggy_beta_', num2str(beta));
suffix_transmittance = strcat('_transmittance_beta_', num2str(beta));

result_format_depth = '.mat';
result_format_foggy = '.png';
result_format_transmittance = '.png';

% Window size for atmospheric light estimation.
window_size = 15;

% ------------------------------------------------------------------------------

% Steps 1 and 2 of the depth denoising and completion pipeline:
% 1. calculation of a raw depth map in meters from the raw disparity input
% 2. denoising and completion of the above raw depth map to produce a refined
%    depth map in meters

% Compute refined depth map.
depth_map =...
    depth_in_meters_cityscapes_stereoscopic_inpainting(left_disparity,...
    camera_parameters_file, left_image, left_image_uint8, right_image);

% Write result to .MAT file.
depth_result_filename = strcat(image_basename, suffix_depth,...
    result_format_depth);
depth_result_dir = fullfile(results_root_dir, depth_results_dir);
if ~exist(depth_result_dir, 'dir')
    mkdir(depth_result_dir);
end
save(fullfile(depth_result_dir, depth_result_filename), 'depth_map');

% ------------------------------------------------------------------------------

% Steps 3 and 4 of the depth denoising and completion pipeline:
% 3. calculation of a scene distance map in meters from the refined depth map
% 4. application of the transmittance formula for a homogeneous medium to obtain
%    an initial transmittance map

t_initial = transmission_homogeneous_medium(depth_map, beta,...
    camera_parameters_file);

% ------------------------------------------------------------------------------

% Step 5 of the depth denoising and completion pipeline:
% 5. guided filtering of the initial transmittance map using the original
%    clear-weather image as guidance to compute the final transmittance map

% Compute final transmittance map.
t = transmission_postprocessing_guided_filter(t_initial, left_image);

% Write result to PNG image.
transmittance_result_filename = strcat(image_basename, suffix_transmittance,...
    result_format_transmittance);
transmittance_result_dir = fullfile(results_root_dir,...
    transmittance_results_dir);
if ~exist(transmittance_result_dir, 'dir')
    mkdir(transmittance_result_dir);
end
imwrite(t, fullfile(transmittance_result_dir, transmittance_result_filename));

% ------------------------------------------------------------------------------

% Estimation of atmospheric light from the clear-weather image, using the method
% proposed by He et al. in "Single Image Haze Removal Using Dark Channel Prior"
% (IEEE T-PAMI, 2011) with the improvement of Tang et al. in "Investigating
% Haze-relevant Features in a Learning Framework for Image Dehazing" (CVPR,
% 2014).

left_image_dark_channel = get_dark_channel(left_image, window_size);
L_atm = estimate_atmospheric_light_rf(left_image_dark_channel, left_image);

% ------------------------------------------------------------------------------

% Fog simulation from the clear-weather image and the estimated transmittance
% map and atmospheric light, using the standard optical model for fog and haze
% introduced by Koschmieder in "Theorie der horizontalen Sichtweite" (Beitrage
% zur Physik der freien Atmosphaere, 1924) which assumes homogeneous atmosphere
% and globally constant atmospheric light.

% Compute partially synthetic foggy image.
I = haze_linear(left_image, t, L_atm);

% Write result to PNG image.
foggy_result_filename = strcat(image_basename, suffix_foggy,...
    result_format_foggy);
foggy_result_dir = fullfile(results_root_dir,foggy_results_dir);
if ~exist(foggy_result_dir, 'dir')
    mkdir(foggy_result_dir);
end
imwrite(I, fullfile(foggy_result_dir, foggy_result_filename));
