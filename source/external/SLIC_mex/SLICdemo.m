%======================================================================
%SLIC demo
% Copyright (C) 2015 Ecole Polytechnique Federale de Lausanne
% File created by Radhakrishna Achanta
% Please also read the copyright notice in the file slicmex.c 
%======================================================================
%Input parameters are:
%[1] 8 bit images (color or grayscale)
%[2] Number of required superpixels (optional, default is 200)
%[3] Compactness factor (optional, default is 10)
%
%Ouputs are:
%[1] labels (in raster scan order)
%[2] number of labels in the image (same as the number of returned
%superpixels
%
%NOTES:
%[1] number of returned superpixels may be different from the input
%number of superpixels.
%[2] you must compile the C file using mex slicmex.c before using the code
%below
%======================================================================
%img = imread('someimage.jpg');
% img = imread('bee.jpg');
% [labels, numlabels] = slicmex(img,500,20);%numlabels is the same as number of superpixels
img = imread('/srv/glusterfs/daid/data/cityscape/leftImg8bit/train/aachen/aachen_000001_000019_leftImg8bit.png');
[labels, numlabels] = slicmex(img,2000,10);%numlabels is the same as number of superpixels
imagesc(labels);
