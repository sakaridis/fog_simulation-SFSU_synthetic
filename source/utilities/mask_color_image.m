function I_masked = mask_color_image(I, binary_mask, mask_color)
%MASK_COLOR_IMAGE  Mask some of the pixels of a color image with a constant
%color.
%
%   INPUTS:
%
%   -|I|: input RGB image.
%
%   -|binary_mask|: logical matrix of the same number of rows and columns as
%   |I|.
%
%   -|mask_color|: 1-by-3 vector containing the RGB color of the applied mask.
%
%   OUTPUTS:
%
%   -|I_masked|: output RGB image, in which pixels where |binary_mask| is true
%   are painted with the |mask_color|.

I_masked_red = I(:, :, 1);
I_masked_red(binary_mask) = mask_color(1);
I_masked_green = I(:, :, 2);
I_masked_green(binary_mask) = mask_color(2);
I_masked_blue = I(:, :, 3);
I_masked_blue(binary_mask) = mask_color(3);
I_masked = cat(3, I_masked_red, I_masked_green, I_masked_blue);

end

