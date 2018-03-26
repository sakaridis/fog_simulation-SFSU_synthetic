function [L, index_L] = estimate_atmospheric_light_dcp(I_dark, I)
%ESTIMATE_ATMOSPHERIC_LIGHT_DCP  Estimate atmospheric light from a fraction of
%brightest pixels in the dark channel of the input image, as proposed in Single
%Image Haze Removal Using Dark Channel Prior.
%
%   INPUTS:
%
%   -|I_dark|: grayscale image of dark channel.
%
%   -|I|: color image of the same height and width as |I_dark|.
%
%   OUTPUTS:
%
%   -|L|: 1-by-1-by-3 matrix with estimated atmospheric light value.
%
%   -|index_L|: linear index for the single-channel version of the image of the
%   pixel with equal value to atmospheric light.

% Determine the number of brightest pixels in dark channel that are used to
% estimate the atmospheric light.
brightest_pixels_fraction = 1 / 1000;
[height, width] = size(I_dark);
number_of_pixels = height * width;
brightest_pixels_count = floor(brightest_pixels_fraction * number_of_pixels);

% Identify indices of brightest pixels in dark channel.
I_dark_vector = reshape(I_dark, number_of_pixels, 1);
[~, indices] = sort(I_dark_vector, 'descend');
brightest_pixels_indices = indices(1:brightest_pixels_count);

% Compute graylevel intensities in original image of brightest pixels in dark
% channel.
I_gray_vector = reshape(rgb2gray(I), number_of_pixels, 1);
I_gray_vector_brightest_pixels = I_gray_vector(brightest_pixels_indices);

% Recover the subscript of the pixel which gives the atmospheric light.
[~, index_highest_intensity] = max(I_gray_vector_brightest_pixels);
index_L = brightest_pixels_indices(index_highest_intensity);
[row_L, column_L] = ind2sub([height, width], index_L);

L = I(row_L, column_L, :);

end

