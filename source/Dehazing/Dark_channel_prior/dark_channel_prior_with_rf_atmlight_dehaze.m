function R = dark_channel_prior_with_rf_atmlight_dehaze(I)
%DARK_CHANNEL_PRIOR_DEHAZE_WITH_RF_ATMLIGHT_DEHAZE  Dehaze input image following
%the pipeline proposed in "Single Image Haze Removal Using Dark Channel Prior",
%with the modification of the paper "Investigating Haze-Relevant Features in a
%Learning Framework for Image Dehazing" for estimating the atmospheric light.
%
%   INPUTS:
%
%   -|I|: H-by-W-by-|image_channels| hazy image.
%
%   OUTPUTS:
%
%   -|R|: output dehazed image representing the estimate for the true radiance
%   of the scene, with same size as |I|.

% Add required paths.
current_script_full_name = mfilename('fullpath');
current_script_directory = fileparts(current_script_full_name);
addpath(fullfile(current_script_directory, '..', '..', 'utilities'));
addpath_relative_to_caller(current_script_full_name,...
    fullfile('..', '..', 'Fog_simulation'));
addpath_relative_to_caller(current_script_full_name,...
    fullfile('..', '..', 'Dehazing'));

% Set parameters.
neighborhood_size_dark_channel = 15;
t_thresh = 0.1;
window_size_guided_filter = 41;
epsilon = 1e-3;

% Dehaze using Dark Channel Prior, following the paper of He, Sun and Tang (CVPR
% 2009). The refinement of the raw transmission map is performed according to
% Guided Image Filtering by He, Sun and Tang (ECCV 2010), in order to gain
% speed.
[I_dark, I_eroded] = get_dark_channel(I, neighborhood_size_dark_channel);
L = estimate_atmospheric_light_rf(I_dark, I);
t_initial = transmission_initial(I_eroded, L);
t = transmission_guided_filtering(t_initial, I,...
    window_size_guided_filter, epsilon);
t = clip_to_unit_range(t);
R = inverse_haze_linear(I, t, L, t_thresh);
R = clip_to_unit_range(R);

end

