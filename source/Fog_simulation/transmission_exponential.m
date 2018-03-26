function t = transmission_exponential(d, beta, varargin)
%TRANSMISSION_EXPONENTIAL  Compute transmission map using given depth map based
%on the Lambert-Beer law
%   Inputs:
%       -|d|: H-by-W matrix with values of depth for processed image in meters.
%       -|beta|: positive constant that parameterizes the density of haze.
%       Larger values indicate thicker haze.
%
%   Outputs:
%       -|t|: H-by-W matrix with medium transmission values ranging in [0, 1].

t = exp(-beta * d);

end

