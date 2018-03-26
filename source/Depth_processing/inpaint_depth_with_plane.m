function inpainted_depth = inpaint_depth_with_plane(pixel_coordinates, plane)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

number_of_pixels = size(pixel_coordinates, 1);

% At each pixel, return the **rectified** interpolation using the plane that is
% assumed to contain all input pixels, in order to always output positive depth
% values.
inpainted_depth = max([pixel_coordinates, ones(number_of_pixels, 1)] * plane,...
    0);

end

