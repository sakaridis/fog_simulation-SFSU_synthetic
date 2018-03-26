function t_refined = transmission_guided_filtering(t, I, window_size, epsilon)
%TRANSMISSION_GUIDED_FILTERING  Applies guided image filter to initial
%transmission map using the input hazy or clean image as guidance image.
%
%   INPUTS:
%
%   -|t|: grayscale image of the same type, height and width as the input image
%   |I|, containing the initial estimate for the transmission map.
%
%   -|I|: H-by-W-by-|image_channels| hazy image.
%
%   -|window_size|: side of square window used by the guided filter, in pixels.
%
%   -|epsilon|: scalar that controls the amount of regularization for the guided
%   filter.
%
%   OUTPUTS:
%
%   -|t_refined|: grayscale image of the same dimensions as |t|, containing the
%   refined transmission map after guided image filtering.

t_refined = imguidedfilter(t, I, 'NeighborhoodSize', window_size,...
    'DegreeOfSmoothing', epsilon);

end

