function disparity_in_pixels = disparity_in_pixels_cityscapes(input_disparity)

disparity_in_pixels = (double(input_disparity) - 1) / 256;

end

