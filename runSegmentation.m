% This script segments biofilms into its single cell constituents. The 
% algorithm takes in a deconvolved image with background removed (and set 
% to zero) and outputs a list of 3d coordinates corresponding to each cell.
%
% Contact Jing Yan, jing.yan@yale.edu, MCDB, Yale University for
% further information.


% Path for the image file which is to be segmented
file_name = 'C:\...\segmentationCode\sample.tif';

% Read the image file
img = imread(file_name,1);
info = imfinfo(file_name);
I = zeros(size(img,1), size(img,2), length(info));
for i=1:length(info)
    img = imread(file_name,i);
    I(:,:,i)=double(img);
end

% Initial segmentation identifying all connected components
[img_binary,img_label,S] = segmentation_cellfilter_YanLab(I);

% Recursive cell segmentation using adaptive thresholding
[S] = segmentation_3Dprocess(img_label,I,S);

