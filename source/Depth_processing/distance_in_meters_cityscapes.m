function distance_map_in_meters = distance_in_meters_cityscapes(...
    depth_map_in_meters, camera_parameters_file)
%DISTANCE_IN_METERS_CITYSCAPES  Compute air thickness, i.e. distance of depicted
%object from camera center, expressed in meters, as a dense map with same
%resolution as the image, using a dense depth map and intrinsic camera
%parameters as inputs.

% Retrieve relevant intrinsic camera parameters: focal length and optical
% center, both expressed in pixel coordinates.
[~, f_x, c_x, c_y] = camera_parameters_cityscapes(camera_parameters_file);

% Compute medium (i.e. air) thickness in meters from depth. The derivation of
% the formula in the final lines of code is based on similar triangles.

[height, width] = size(depth_map_in_meters);
[X, Y] = meshgrid(1:width, 1:height);

distance_map_in_meters = depth_map_in_meters .*...
    sqrt((f_x ^ 2 + (X - c_x) .^ 2 + (Y - c_y) .^ 2) / f_x ^ 2);

end

