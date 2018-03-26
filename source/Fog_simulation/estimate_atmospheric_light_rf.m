function [L, index_L] = estimate_atmospheric_light_rf(I_dark, I)
%ESTIMATE_ATMOSPHERIC_LIGHT_RF  Estimate atmospheric light from a fraction of
%brightest pixels in the dark channel of the input image, as proposed in
%Investigating Haze-Relevant Features in a Learning Framework for Image
%Dehazing.
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
% estimate the atmospheric light. Ensure that this number is always odd, so that
% the final median filter always outputs a value that is equal to the intensity
% of a pixel.
brightest_pixels_fraction = 1 / 1000;
[height, width] = size(I_dark);
number_of_pixels = height * width;
brightest_pixels_count = brightest_pixels_count_rf(number_of_pixels,...
    brightest_pixels_fraction);

% Identify indices of brightest pixels in dark channel.
I_dark_vector = reshape(I_dark, number_of_pixels, 1);
[~, indices] = sort(I_dark_vector, 'descend');
brightest_pixels_indices = indices(1:brightest_pixels_count);

% Compute graylevel intensities in original image of brightest pixels in dark
% channel.
I_gray_vector = reshape(rgb2gray(I), number_of_pixels, 1);
I_gray_vector_brightest_pixels = I_gray_vector(brightest_pixels_indices);

% Recover the subscript of the pixel which gives the atmospheric light as the
% one from the brightest pixels with median graylevel intensity in original
% image.
median_intensity = median(I_gray_vector_brightest_pixels);
index_median_intensity =...
    find(I_gray_vector_brightest_pixels == median_intensity, 1);
index_L = brightest_pixels_indices(index_median_intensity);
[row_L, column_L] = ind2sub([height, width], index_L);

L = I(row_L, column_L, :);

end

