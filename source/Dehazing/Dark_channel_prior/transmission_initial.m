function t = transmission_initial(I_eroded, L)
%TRANSMISSION_INITIAL  Compute transmission map using dark channel of normalized
%hazy image.
%
%   INPUTS:
%
%   -|I_eroded|: intermediate eroded image of the same dimensions as |I|. It is
%   used to avoid recomputations, as the estimation of atmospheric light has
%   already required an erosion of the input image.
%
%   -|L|: 1-by-1-by-3 matrix with estimated atmospheric light value.
%
%   OUTPUTS:
%
%   -|t|: grayscale image of the same type, height and width as the input image
%   |I|, containing the initial estimate for the transmission map.

% Parameter controlling the amount of haze that is removed.
omega = 0.95;

t = 1 - omega * min(I_eroded ./...
    repmat(L, size(I_eroded, 1), size(I_eroded, 2)), [], 3);

end

