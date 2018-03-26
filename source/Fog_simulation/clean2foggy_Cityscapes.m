function clean2foggy_Cityscapes(left_image_file_names,...
    right_image_file_names, disparity_file_names,...
    camera_parameters_file_names, depth_file_names, load_depth, save_depth,...
    depth_completion_method, attenuation_coefficient_method, beta_parameters,...
    atmospheric_light_estimation, transmittance_model,...
    transmittance_postprocessing, fog_optical_model,...
    cityscapes_foggy_output_directory,...
    cityscapes_transmittance_output_directory, output_format)
%CLEAN2FOGGY_CITYSCAPES  Create foggy counterparts for an input list of clean
%Cityscapes images.
%
%   INPUTS:
%
%   -|left_image_file_names|: cell array of strings with full paths to clean
%    left input images.
%
%   -|right_image_file_names|: cell array of strings with full paths to clean
%    right input images.
%
%   -|disparity_file_names|: cell array of strings with full paths to disparity
%    maps for input images.
%
%   -|camera_parameters_file_names|: cell array of strings with full paths to
%    JSON files where camera parameters are stored.
%
%   -|depth_file_names|: cell array of strings with full paths to depth maps for
%    input images, whether these depth map files exist or not.
%
%   -|load_depth|: flag indicating whether depth should be loaded from existing
%    depth map files if possible.
%
%   -|save_depth|: flag indicating whether the computed depth maps should be
%    saved to files to accelerate future runs.
%
%   -|depth_completion_method|: handle of the function which implements
%    denoising and completion of the raw depth map.
%
%   -|attenuation_coefficient_method|: handle of the function which is specified
%    as the attenuation coefficient method.
%
%   -|beta_parameters|: struct with fields that help compute values for
%    attenuation coefficient beta corresponding to the processed images.
%
%   -|atmospheric_light_estimation|: handle of the function which is used to
%    estimate atmospheric light.
%
%   -|transmittance_model|: handle of the function that implements computation
%    of transmittance map given a depth map.
%
%   -|transmittance_postprocessing|: handle of the function which is used to
%    postprocess the transmittance map after the first denoising and completion
%    step.
%
%   -|fog_optical_model|: handle of the function that implements fog simulation
%    on a clean image.
%
%   -|cityscapes_foggy_output_directory|: full path to root directory where the
%    synthetic foggy output Cityscapes images are saved.
%
%   -|cityscapes_transmittance_output_directory|: full path to root directory
%    where the transmittance maps for foggy Cityscapes images are saved.
%
%   -|output_format|: string specifying the image format for the foggy output
%    images, e.g. '.png'

% Total number of processed images. Should be equal to number of files for each
% auxiliary set.
number_of_images = length(left_image_file_names);
assert(number_of_images == length(right_image_file_names));
assert(number_of_images == length(disparity_file_names));
assert(number_of_images == length(camera_parameters_file_names));
assert(number_of_images == length(depth_file_names));

% Generate attenuation coefficient values for the set of processed images.
beta_vector = attenuation_coefficient_method(number_of_images, beta_parameters);

% Compute and save foggy images and transmittance maps. Save intermediately
% computed depth map optionally.
for i = 1:number_of_images
    
    % Read left image and bring it to double precision for subsequent
    % computations.
    R_left_uint8 = imread(left_image_file_names{i});
    R_left = im2double(R_left_uint8);
    
    % Check whether fog simulation has already been run for current image at an
    % earlier point in time.
    result_exists = exist(depth_file_names{i}, 'file');
    
    % Get depth. Store it in |depth_map|, which is a 2-dimensional matrix in
    % double format with the same resolution as the input image, containing
    % depth values in meters.
    load_depth_current = load_depth && result_exists;
    if load_depth_current
        % Load depth into variable |depth_map|.
        load(depth_file_names{i}, 'depth_map');
        
    else
        % Use both images of the stereo pair, the disparity map for the left
        % view and the intrinsic parameters of the camera to compute a complete
        % depth map (possibly denoised).
        
        % Read input disparity map.
        input_disparity = imread(disparity_file_names{i});
        
        % Read right image of stereo pair to double precision.
        R_right = im2double(imread(right_image_file_names{i}));
    
        % Compute depth.
        depth_map = depth_completion_method(input_disparity,...
            camera_parameters_file_names{i}, R_left, R_left_uint8, R_right);
    end
    
    % Prevent saving the depth map when the corresponding .mat file already
    % exists and is being loaded.
    save_depth_current = save_depth && ~result_exists;
    
    % If the depth map should be saved, do it after creating the appropriate
    % directory.
    if save_depth_current
        current_depth_output_directory = fileparts(depth_file_names{i});
        if exist(current_depth_output_directory) ~= 7
            mkdir(current_depth_output_directory);
        end
        save(depth_file_names{i}, 'depth_map');
    end

    % Compute transmittance map using the specified transmittance model.
    t = transmittance_model(depth_map, beta_vector(i),...
        camera_parameters_file_names{i});
    
    % Postprocessing of transmittance map. May stand for no postprocessing.
    t = transmittance_postprocessing(t, R_left);
    
    % Estimate atmospheric light for the clean image.
    neighborhood_size_dark_channel = 15;
    R_left_dark = get_dark_channel(R_left, neighborhood_size_dark_channel);
    L_atm = atmospheric_light_estimation(R_left_dark, R_left);
    
    % Simulate fog using the specified fog optical model and the computed
    % transmittance map and atmospheric light for the current image.
    I = fog_optical_model(R_left, t, L_atm);
    
    % --------------------------------------------------------------------------
    
    % Save the foggy image and the transmittance map in the respective output
    % directories in a (preferably) lossless format.
    
    % The name of and path to the output image is based on the input image.
    [path_to_input, R_left_name] =...
        fileparts(left_image_file_names{i});
    
    % Suffices specifying the parameter values that were used to generate haze.
    parameters_suffix_foggy = strcat('_foggy_beta_', num2str(beta_vector(i)));
    parameters_suffix_transmittance = strcat('_transmittance_beta_',...
        num2str(beta_vector(i)));
    
    % Full names of output images formed by the name of the input image and the
    % suffices.
    I_name_with_extension = strcat(R_left_name, parameters_suffix_foggy,...
        output_format);
    t_name_with_extension = strcat(R_left_name,...
        parameters_suffix_transmittance, output_format);
    
    % Determine output directories based on the directory structure of original
    % Cityscapes dataset: 1) train-val-test directories, 2) city directories.
    path_to_input_split = strsplit(path_to_input, filesep);
    current_foggy_output_directory =...
        fullfile(cityscapes_foggy_output_directory,...
        path_to_input_split{end - 1}, path_to_input_split{end});
    current_transmittance_output_directory =...
        fullfile(cityscapes_transmittance_output_directory,...
        path_to_input_split{end - 1}, path_to_input_split{end});
    
    % Create output directories where synthetic foggy images and transmittance
    % maps will be saved, if they do not already exist.
    if exist(current_foggy_output_directory) ~= 7
        mkdir(current_foggy_output_directory);
    end
    if exist(current_transmittance_output_directory) ~= 7
        mkdir(current_transmittance_output_directory);
    end
    
    % Build image names from base name, parameters suffices and output format
    % and save them if not already saved.
    I_exists = exist(fullfile(current_foggy_output_directory,...
        I_name_with_extension), 'file');
    if ~I_exists
        imwrite(I,...
            fullfile(current_foggy_output_directory, I_name_with_extension));
    end
    
    t_exists = exist(fullfile(current_transmittance_output_directory,...
        t_name_with_extension), 'file');
    if ~t_exists
        imwrite(t, fullfile(current_transmittance_output_directory,...
            t_name_with_extension));
    end
    
end

end

