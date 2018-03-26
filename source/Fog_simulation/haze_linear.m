function I = haze_linear(R, t, L)
%HAZE_LINEAR  Generate hazy image from clean image using the linear haze model
%corresponding to Lambert-Beer law
%   Inputs:
%       -|R|: H-by-W-by-|image_channels| clean image representing true radiance
%       of scene.
%       -|t|: H-by-W transmission map.
%       -|L|: 1-by-1-by-|image_channels| homogeneous atmospheric light.
%
%   Outputs:
%       -|I|: synthetic hazy image, with same size as input clean image |R|.

image_channels = size(L, 3);

% Auxiliary matrix with replicates of transmission map for all channels,
% allowing to express hazy image conveniently.
t_replicated = repmat(t, 1, 1, image_channels);

% Apply linear haze model.
I = t_replicated .* R + (1 - t_replicated) .* repmat(L, size(t));

end

