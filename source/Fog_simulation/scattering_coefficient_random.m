function beta_vector = scattering_coefficient_random(number_of_images,...
    parameters)
%SCATTERING_COEFFICIENT_RANDOM  Generate scattering coefficient values for a set
%of images uniformly at random.
%   Inputs:
%       -|number_of_images|: number of images for which scattering coefficient
%       is simulated.
%       -|parameters|: structure containing miscellaneous parameters, such as
%       range of values for scattering coefficient and type of random number
%       generator. Guarantees uniformity with other functions that implement
%       generation of scattering coefficient values.
%
%   Outputs:
%       -|beta_vector|: 1-by-|number_of_images| vector containing the random
%       values of scattering coefficient for every image in the set.

% Determine the range of random values for scattering coefficient.
maximum_value = parameters.maximum_value;
minimum_value = parameters.minimum_value;

% Get type of random number generator, e.g. 'default'.
random_generator = parameters.random_generator;

% Get binary flag that indicates whether the random number generator should be
% configured or not.
configure_random_generator = parameters.configure_random_generator;

% Optionally configure random number generation for repeatability.
if configure_random_generator
    rng(random_generator);
end

% Generate random scattering coefficient for each image following a uniform
% distribution inside the specified range.
beta_vector = minimum_value + (maximum_value - minimum_value) *...
    rand(1, number_of_images);

end

