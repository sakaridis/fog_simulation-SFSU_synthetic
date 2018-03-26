function beta_vector = scattering_coefficient_fixed(number_of_images,...
    parameters)
%SCATTERING_COEFFICIENT_FIXED  Generate fixed scattering coefficient for a set
%of images.
%   Inputs:
%       -|number_of_images|: number of images for which scattering coefficient
%       is simulated.
%       -|parameters|: structure containing miscellaneous parameters, such as
%       fixed scattering coefficient value. Guarantees uniformity with other
%       functions that implement generation of scattering coefficient values.
%
%   Outputs:
%       -|beta_vector|: 1-by-|number_of_images| vector containing the random
%       values of scattering coefficient for every image in the set.

% Determine the fixed intensity of atmospheric light.
beta = parameters.beta;

% Generate fixed scattering coefficient for all images.
beta_vector = repmat(beta, 1, number_of_images);

end

