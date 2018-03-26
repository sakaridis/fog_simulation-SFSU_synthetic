function clipped_image = clip_to_unit_range(image)

clipped_image = min(max(image, 0), 1);

end