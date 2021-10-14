function [img_bac_sharpened_binary_filtered,img_bac_sharpened_label_filtered,S] = segmentation_cellfilter_YanLab(img_bac_sharpened,cell_size_lower_thresh)
% Takes in a deconvolved biofilm z-stack image with background removed and
% identifies all connected components. 
%
% inputs: 
% img_bac_sharpened - deconvolved biofilm with background removed and set
% to 0.
% cell_size_lower_thresh - minimum voxel size of a single bacterium.
% Connected components smaller than this value are removed. 

% output:
% img_bac_sharpened_binary_filtered - binarized image
% img_bac_sharpened_label_filtered - image labelled by connected component
% number
% S - list of all connected components. 

if nargin == 1
    cell_size_lower_thresh = 25;
end

%Binarize the image. 
img_bac_sharpened_binary_unfiltered=img_bac_sharpened;
img_bac_sharpened_binary_unfiltered(img_bac_sharpened_binary_unfiltered(:,:,:)>0)=1;

% Identify all connected components and populate
% img_bac_sharpened_label_filtered and S, do not populate small cells.  
X=bwconncomp(img_bac_sharpened_binary_unfiltered(:,:,1:size(img_bac_sharpened_binary_unfiltered,3)),6); 
XX=X.PixelIdxList;
num_of_clusters=1;
xlength=size(img_bac_sharpened,1);
ylength=size(img_bac_sharpened,2);
S=([]);
img_bac_sharpened_label_filtered=0*img_bac_sharpened;
for nofclustersfake=1:size(XX,2)
    if size(XX{nofclustersfake},1)>cell_size_lower_thresh
        S{num_of_clusters}=zeros(size(XX{nofclustersfake},1)+1,3);
        S{num_of_clusters}(1,1)=size(XX{nofclustersfake},1);
        for i=1:size(XX{nofclustersfake},1)
            S{num_of_clusters}(i+1,3)=1+floor(XX{nofclustersfake}(i,1)/(xlength*ylength+0.00001));
            S{num_of_clusters}(i+1,2)=1+floor((XX{nofclustersfake}(i,1)-(S{num_of_clusters}(i+1,3)-1)*(xlength*ylength))/(xlength+0.00001));
            S{num_of_clusters}(i+1,1)=XX{nofclustersfake}(i,1)-(S{num_of_clusters}(i+1,3)-1)*(xlength*ylength)-(S{num_of_clusters}(i+1,2)-1)*(xlength);
            img_bac_sharpened_label_filtered(S{num_of_clusters}(i+1,1),S{num_of_clusters}(i+1,2),S{num_of_clusters}(i+1,3))=num_of_clusters;
        end
        
        num_of_clusters=num_of_clusters+1;
    end
end

% Create the binarized image without small cells. 
img_bac_sharpened_binary_filtered=img_bac_sharpened_label_filtered;
img_bac_sharpened_binary_filtered(img_bac_sharpened_binary_filtered(:,:,:)>0)=1;

num_of_clusters = num_of_clusters-1;% output total number of cells as a quick check.
disp(['number of clusters :',num2str(num_of_clusters)])




