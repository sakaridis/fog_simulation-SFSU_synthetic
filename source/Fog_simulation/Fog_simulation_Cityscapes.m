function Fog_simulation_Cityscapes(task_id, dataset_split, refinement_level,...
    variant, beta, output_root_directory, images_per_task)
%FOG_SIMULATION_CITYSCAPES  Simulate fog for a batch of images from Cityscapes
%and write the results. Structured for execution on a cluster.
%
%   INPUTS:
%
%   -|task_id|: ID of the task. Used in conjuction with |images_per_task| to
%    determine which images out of the entire dataset will form the batch that
%    will be processed by this task.
%
%   -|dataset_split|: string that indicates which subset of Cityscapes is used
%    for running the simulation, e.g. 'trainval', 'train_extra'.
%
%   -|refinement_level|: string that indicates the refinement that is applied to
%    the used subset, e.g. 'full', 'refined'.
%
%   -|variant|: string that indicates the type of fog simulation that is being
%    run, e.g. 'stereoscopic_inpainting_with_guided_filtering',
%    'nearest_neighbor'.
%
%   -|beta|: value of attenuation coefficient that determines fog density.
%
%   -|cityscapes_output_root_directory|: full path to directory under which the
%    results of the simulation are written, e.g.
%    '/scratch_net/nowin/csakarid/Toyota-foggy/data/Cityscapes'.
%
%   -|images_per_task|: maximum number of images for which fog simulation is run
%    in each task.
%    19997 train_extra Cityscapes images -> 200 images per task * 100 tasks.
%    5000 trainvaltest Cityscapes images -> 100 images per task * 50 tasks.
%    550 trainval refined Cityscapes images -> 22 images per task * 25 tasks.

if ischar(task_id)
    task_id = str2double(task_id);
end

% ------------------------------------------------------------------------------

% Check which variant of the simulation is being run and distinguish between
% different cases.

% Get abbreviations to be included in output names.
abbreviations_of_variants =...
    containers.Map({'stereoscopic_inpainting_with_guided_filtering',...
    'nearest_neighbor'}, {'stereogf', 'nearest'});
variant_abbreviation = abbreviations_of_variants(variant);

% ------------------------------------------------------------------------------

% Add paths to functions that are called for fog simulation.

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

% ------------------------------------------------------------------------------

% Create lists of input files for various modalities.

cityscapes_file_lists_directory = fullfile(current_script_directory, '..',...
    '..', 'data', 'lists_file_names');

% Left input images.
list_of_left_images_file = strcat('leftImg8bit_', dataset_split, '_',...
    refinement_level, '_clean_filenames.txt');
fid = fopen(fullfile(cityscapes_file_lists_directory,...
    list_of_left_images_file));
left_image_file_names = textscan(fid, '%s');
fclose(fid);
left_image_file_names = left_image_file_names{1};

number_of_images = length(left_image_file_names);

% Right input images.
list_of_right_images_file = strcat('rightImg8bit_', dataset_split, '_',...
    refinement_level, '_filenames.txt');
fid = fopen(fullfile(cityscapes_file_lists_directory,...
    list_of_right_images_file));
right_image_file_names = textscan(fid, '%s');
fclose(fid);
right_image_file_names = right_image_file_names{1};

% Disparity map files.
list_of_disparity_maps_file = strcat('disparity_', dataset_split, '_',...
    refinement_level, '_filenames.txt');
fid = fopen(fullfile(cityscapes_file_lists_directory,...
    list_of_disparity_maps_file));
disparity_file_names = textscan(fid, '%s');
fclose(fid);
disparity_file_names = disparity_file_names{1};

% Camera parameters files.
list_of_camera_parameters_file = strcat('camera_', dataset_split, '_',...
    refinement_level, '_filenames.txt');
fid = fopen(fullfile(cityscapes_file_lists_directory,...
    list_of_camera_parameters_file));
camera_parameters_file_names = textscan(fid, '%s');
fclose(fid);
camera_parameters_file_names = camera_parameters_file_names{1};

% Based on the selected simulation variant, form the name of the corresponding
% text file with the file names of the respective depth maps.
switch variant_abbreviation
    case 'stereogf'
        list_of_depth_maps_file = strcat('depth_stereoscopic_',...
            dataset_split, '_', refinement_level, '_filenames.txt');
    case 'nearest'
        list_of_depth_maps_file = strcat('depth_nearest_', dataset_split,...
            '_', refinement_level, '_filenames.txt');
end

% Get all names of depth map files, whether they already exist or not.
fid = fopen(fullfile(cityscapes_file_lists_directory,...
    list_of_depth_maps_file));
depth_file_names = textscan(fid, '%s');
fclose(fid);
depth_file_names = depth_file_names{1};

% Sanity checks. Total number of files should be identical for all modalities.
assert(number_of_images == length(right_image_file_names));
assert(number_of_images == length(disparity_file_names));
assert(number_of_images == length(camera_parameters_file_names));
assert(number_of_images == length(depth_file_names));

% ------------------------------------------------------------------------------

% Determine current batch.

% Determine the set of images that are to be processed in the current task.
batch_ind = (task_id - 1) * images_per_task + 1:task_id * images_per_task; 
if batch_ind(1) > number_of_images
    return;
end
if batch_ind(end) > number_of_images
    % Truncate for last task.
    batch_ind = batch_ind(1:number_of_images - batch_ind(1) + 1);
end

% ------------------------------------------------------------------------------

% Instantiate individual components of the fog simulation pipeline that are
% required for the selected variant and set their parameters.

% Attenuation coefficient. Set it to meaningful values.
attenuation_coefficient_method = @scattering_coefficient_fixed;
beta_parameters.beta = beta;

% Atmospheric light.
atmospheric_light_estimation = @estimate_atmospheric_light_rf;

% Select the optical model which is used for the synthesis of foggy images from
% the clean input images.
fog_optical_model = @haze_linear;

% Transmittance model of the clear scene radiance.
transmittance_model = @transmission_homogeneous_medium;

% Specify whether depth maps will be loaded from existing .mat files if possible
% or computed on the fly.
load_depth = true;

% Specify whether the depth maps used for fog simulation will be saved to .mat
% files.
save_depth = true;

% Select the method for depth map completion and denoising, and the method for
% subsequent postprocessing of the initial transmittance map.
switch variant_abbreviation
    case 'stereogf'
        depth_completion_method =...
            @depth_in_meters_cityscapes_stereoscopic_inpainting;
        transmittance_postprocessing =...
            @transmission_postprocessing_guided_filter;
    case 'nearest'
        depth_completion_method = @depth_in_meters_cityscapes_nn_inpainting;
end

% ------------------------------------------------------------------------------

% Output specifications.

% Specify output format for synthetic foggy images and corresponding
% transmittance maps.
output_format = '.png';

% Determine initial part of the path of output files.
cityscapes_output_directories_basename = strcat('leftImg8bit_',...
    dataset_split, '_', refinement_level, '_', variant_abbreviation,...
    '_beta_', num2str(beta_parameters.beta));
cityscapes_output_foggy_directory = fullfile(output_root_directory,...
    strcat(cityscapes_output_directories_basename, '_foggy'));
cityscapes_output_transmittance_directory = fullfile(output_root_directory,...
    strcat(cityscapes_output_directories_basename, '_transmittance'));

% ------------------------------------------------------------------------------

% Run fog simulation on current batch using the specified settings.

switch variant_abbreviation
    
    % Stereoscopic inpainting followed by guided filtering.
    case 'stereogf'
        clean2foggy_Cityscapes(left_image_file_names(batch_ind),...
            right_image_file_names(batch_ind),...
            disparity_file_names(batch_ind),...
            camera_parameters_file_names(batch_ind),...
            depth_file_names(batch_ind), load_depth, save_depth,...
            depth_completion_method, attenuation_coefficient_method,...
            beta_parameters, atmospheric_light_estimation,...
            transmittance_model, transmittance_postprocessing,...
            fog_optical_model, output_root_directory,...
            cityscapes_output_foggy_directory,...
            cityscapes_output_transmittance_directory, output_format);
        
    % Nearest neighbor inpainting of depth values.
    case 'nearest'
        clean2foggy_Cityscapes_nn(left_image_file_names(batch_ind),...
            disparity_file_names(batch_ind),...
            camera_parameters_file_names(batch_ind),...
            attenuation_coefficient_method, beta_parameters,...
            atmospheric_light_estimation, transmittance_model,...
            fog_optical_model, cityscapes_output_foggy_directory,...
            cityscapes_output_transmittance_directory, output_format);
end

end

