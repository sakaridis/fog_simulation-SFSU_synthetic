function [I_dark, I_eroded] = get_dark_channel(I, neighborhood_size)
%GET_DARK_CHANNEL  Compute dark channel of input image with respect to a square
%neighborhood patch using erosion.
%
%   INPUTS:
%
%   -|I|: input color or grayscale image.
%
%   -|neighborhood_size|: the side of the square patch used for erosion, in
%   pixels.
%
%   OUTPUTS:
%
%   -|I_dark|: output grayscale image of the same type, height and width as |I|.
%
%   -|I_eroded|: intermediate eroded image of the same dimensions as |I|.

se = strel('square', neighborhood_size);

I_eroded = imerode(I, se);
I_dark = min(I_eroded, [], 3);

end