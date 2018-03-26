function image_inpainted = nearest_neighbor_inpainting(image, mask)
%NEAREST_NEIGHBOR_INPAINTING  Inpaint parts of an input image, using
%nearest neighbor inpainting, i.e. copying the values of the respective nearest
%pixels outside the inpainted parts.

% Find the outer boundary of the mask: only pixels in that set are candidates
% for nearest neighbors of pixels in the mask.
mask_dilated = imdilate(mask, strel('square', 3));
mask_outer_boundary = mask_dilated & ~mask;
inds_boundary = find(mask_outer_boundary);

% Get the coordinates of "outer boundary" and "mask" pixels for nearest neighbor
% search.
[height, width] = size(mask);
[X, Y] = meshgrid(1:width, 1:height);
pixel_coords = [X(:), Y(:)];
pixel_coords_boundary = pixel_coords(inds_boundary, :);
pixel_coords_mask = pixel_coords(mask, :);

% Find nearest neighbors for all "mask" pixels.
inds_nearest_neighbors = knnsearch(pixel_coords_boundary, pixel_coords_mask);

% Inpaint pixels in the mask using the values at the nearest neighbor pixels.
image_inpainted = image;
image_inpainted(mask) = image_inpainted(inds_boundary(inds_nearest_neighbors));

end

