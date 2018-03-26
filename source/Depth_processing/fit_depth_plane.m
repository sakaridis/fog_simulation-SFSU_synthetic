function plane = fit_depth_plane(x_hom, depth)
%FIT_DEPTH_PLANE  Least squares fitting of a plane in the 3D scene projected on
%the image from a set of input points.
%
%   INPUTS:
%
%   -|x_hom|: |P|-by-3 matrix holding homogeneous image-plane coordinates of
%   input points.
%
%   -|depth|: |P|-by-1 vector holding depth values of input points.
%
%   OUTPUTS:
%
%   -|plane|: 3-by-1 vector holding the parameters of the fitted plane, so that
%   for each point [x, y, d] it holds d = [x, y, 1] * plane.

x_hom_T = x_hom.';
plane = (x_hom_T * x_hom) \ (x_hom_T * depth);

end

