function R = inverse_haze_linear(I, t, L, t_thresh)
%INVERSE_HAZE_LINEAR  Generate clean image from hazy image using the inversion
%of the linear haze model corresponding to Lambert-Beer law.
%
%   INPUTS:
%
%   -|I|: H-by-W-by-|image_channels| hazy image.
%
%   -|t|: H-by-W transmission map.
%
%   -|L|: 1-by-1-by-|image_channels| homogeneous atmospheric light.
%
%   -|t_thresh|: lower threshold for transmission map values, useful to avoid
%   boosting noise in areas with very low transmission.
%
%   OUTPUTS:
%
%   -|R|: output dehazed image representing the estimate for the true radiance
%   of the scene, with same size as input hazy image |I|.

image_channels = size(L, 3);

% Auxiliary matrix, helping to express clean scene radiance conveniently.
L_replicated = repmat(L, size(t));

% Auxiliary matrix with replicates of transmission map for all channels,
% allowing to express clean scene radiance conveniently.
t_replicated = repmat(t, 1, 1, image_channels);

% Solve equation of homogeneous linear haze model for scene radiance.
R = (I - L_replicated) ./ max(t_replicated, t_thresh) + L_replicated;

end

