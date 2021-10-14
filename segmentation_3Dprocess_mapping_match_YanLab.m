function [match] = segmentation_3Dprocess_mapping_match_YanLab(x,y,z,nofnewcell,onecell_mac_deleterepeat,isearch)
% This function is a sub-program of segmentation_3Dprocess_mapping and is non-executable
% This function determines the belonging of the deleted voxels.
% Input: Information of a deleted voxel
% Output: The belonging of this deleted voxel

distance_min=size(onecell_mac_deleterepeat,1)^2+size(onecell_mac_deleterepeat,2)^2+size(onecell_mac_deleterepeat,3)^2;% Initial amount
nearby_Repaired_cells=[];num=1;

searchrange=1;
while size(nearby_Repaired_cells,1)==0
    % Control the lower search range to avoid border
    if x(isearch)-searchrange>0 
        x_range_1=(x(isearch)-searchrange);
    else
        x_range_1=1;
    end
    if y(isearch)-searchrange>0 
        y_range_1=(y(isearch)-searchrange);
    else
        y_range_1=1;
    end
    if z(isearch)-searchrange>0 
        z_range_1=(z(isearch)-searchrange);
    else
        z_range_1=1;
    end
    % Control the upper search range to avoid border
    if x(isearch)+searchrange<size(onecell_mac_deleterepeat,1)
        x_range_2=(x(isearch)+searchrange);
    else
        x_range_2=size(onecell_mac_deleterepeat,1);
    end
    if y(isearch)+searchrange<size(onecell_mac_deleterepeat,2) 
        y_range_2=(y(isearch)+searchrange);
    else
        y_range_2=size(onecell_mac_deleterepeat,2);
    end
    if z(isearch)+searchrange<size(onecell_mac_deleterepeat,3) 
        z_range_2=(z(isearch)+searchrange);
    else
        z_range_2=size(onecell_mac_deleterepeat,3);
    end
    
    % Continueously enlarge the searching area. This algorithm is not perfect but also relatively low-computation
    for isearch_x=x_range_1:x_range_2
        for isearch_y=y_range_1:y_range_2
            for isearch_z=z_range_1:z_range_2
                if onecell_mac_deleterepeat(isearch_x,isearch_y,isearch_z)>0.6 % All the deleted voxels are marked with 0.5
                    distance=abs(isearch_x-x(isearch))^2+abs(isearch_y-y(isearch))^2+abs(isearch_z-z(isearch))^2; % Calculate the distance between two positions
                    if distance<=distance_min
                        % If equal, jump to next line. If small, clear the old matrix and jump to next line
                        if distance<distance_min
                            distance_min=distance;
                            num=1;nearby_Repaired_cells=[];
                        end
                        if distance==distance_min
                            if num>1 % Two nearby "Repaired cells" share the same Euclid distance
                                if onecell_mac_deleterepeat(isearch_x,isearch_y,isearch_z)~=nearby_Repaired_cells(num-1,1)
                                    nearby_Repaired_cells(num,1)=onecell_mac_deleterepeat(isearch_x,isearch_y,isearch_z);
                                    num=num+1;
                                end
                            else % This repaired cell has the closest Euclid distance
                                nearby_Repaired_cells(num,1)=onecell_mac_deleterepeat(isearch_x,isearch_y,isearch_z);
                                num=num+1;
                            end
                        end
                        
                    end
                end
            end
        end
    end
    searchrange=searchrange+1;
end

match=zeros(num-1,1);% Match is for mark which 'Repaired cells' this deleted voxel should be put into
for isearch_num=1:(num-1)
    match(isearch_num,1)=nearby_Repaired_cells(isearch_num,1);
end