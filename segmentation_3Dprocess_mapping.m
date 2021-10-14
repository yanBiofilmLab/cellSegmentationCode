function [Sreplace,Sadd] = segmentation_3Dprocess_mapping(SS,Sreplace,Sadd,nofnewcell)
% This function is a sub-program of segmentation_3Dprocess and is non-executable
% This function will map each voxel in a cluster that aren't catogorized into any core into a nearby core.
% Input: SS (cluster) Sreplace & Sadd(core)
% Output: Repaired cells 


if nofnewcell>0% Only loop for clusters that split into more than one core 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% Function part 1: Labeling deleted voxels %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Pick out voxels in a cluster that are not in any core
    % Find their coordinates and save them in x, y, z
    xmin=min(SS(2:(SS(1,1)+1),1));xmax=max(SS(2:(SS(1,1)+1),1));
    ymin=min(SS(2:(SS(1,1)+1),2));ymax=max(SS(2:(SS(1,1)+1),2));
    zmin=min(SS(2:(SS(1,1)+1),3));zmax=max(SS(2:(SS(1,1)+1),3));
    
    onecell_mac_deleterepeat=zeros(xmax-xmin+1,ymax-ymin+1,zmax-zmin+1);
    for ilabel_SS=2:(SS(1,1)+1)
        onecell_mac_deleterepeat(SS(ilabel_SS,1)-xmin+1,SS(ilabel_SS,2)-ymin+1,SS(ilabel_SS,3)-zmin+1)=0.5;
    end
    for isearch_Sreplace=2:(Sreplace(1,1)+1)
        onecell_mac_deleterepeat(Sreplace(isearch_Sreplace,1)-xmin+1,Sreplace(isearch_Sreplace,2)-ymin+1,Sreplace(isearch_Sreplace,3)-zmin+1)=1;
    end
    for isearch_newcell=1:nofnewcell
        for isearch_Sadd=2:(Sadd{isearch_newcell}(1,1)+1)
            onecell_mac_deleterepeat(Sadd{isearch_newcell}(isearch_Sadd,1)-xmin+1,Sadd{isearch_newcell}(isearch_Sadd,2)-ymin+1,Sadd{isearch_newcell}(isearch_Sadd,3)-zmin+1)=isearch_newcell+1;
        end
    end
    [x,y,z] = ind2sub(size(onecell_mac_deleterepeat),find(onecell_mac_deleterepeat == 0.5));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% Function part 2: Generating matrix match %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Sadd2=Sadd;Sreplace2=Sreplace;isearch_mum=size(x,1);
    
    % The following method is for shrink the size of the generated matrix. 
    search_area=9000;% Number of processed uncatogorized voxels in one loop 
    search_time=floor(isearch_mum/search_area);% Number of loops to finish the repair
    search_residue=isearch_mum-(search_time)*search_area; 
    
    for isearch_time=1:search_time
        match_middle=cell(1,search_area);
        for isearch=1:search_area
            [match] = segmentation_3Dprocess_mapping_match_YanLab(x,y,z,nofnewcell,onecell_mac_deleterepeat,isearch+search_area*(isearch_time-1));
            match_middle{isearch}=match;
        end
        for isearch=1:search_area
            for isearch_match=1:size(match_middle{isearch},1)% For each number in one voxel's matching process
                matching_number=match_middle{isearch}(isearch_match,1);
                if matching_number==1
                    Sreplace2(1,1)=Sreplace2(1,1)+1;
                    Sreplace2(Sreplace2(1,1)+1,1)=x(isearch+search_area*(isearch_time-1))+xmin-1;
                    Sreplace2(Sreplace2(1,1)+1,2)=y(isearch+search_area*(isearch_time-1))+ymin-1;
                    Sreplace2(Sreplace2(1,1)+1,3)=z(isearch+search_area*(isearch_time-1))+zmin-1;
                end
                
                if matching_number>1
                    Sadd2{matching_number-1}(1,1)=Sadd2{matching_number-1}(1,1)+1;
                    Sadd2{matching_number-1}(Sadd2{matching_number-1}(1,1)+1,1)=x(isearch+search_area*(isearch_time-1))+xmin-1;
                    Sadd2{matching_number-1}(Sadd2{matching_number-1}(1,1)+1,2)=y(isearch+search_area*(isearch_time-1))+ymin-1;
                    Sadd2{matching_number-1}(Sadd2{matching_number-1}(1,1)+1,3)=z(isearch+search_area*(isearch_time-1))+zmin-1;
                end
                
            end
        end
    end
    
    for isearch_time=search_time+1
        match_middle=cell(1,search_residue);
        for isearch=1:search_residue
            [match] = segmentation_3Dprocess_mapping_match_YanLab(x,y,z,nofnewcell,onecell_mac_deleterepeat,isearch+search_area*(isearch_time-1));
            match_middle{isearch}=match;
        end
        for isearch=1:search_residue
            for isearch_match=1:size(match_middle{isearch},1)% For each number in one voxel's match
                matching_number=match_middle{isearch}(isearch_match,1);
                if matching_number==1
                    Sreplace2(1,1)=Sreplace2(1,1)+1;
                    Sreplace2(Sreplace2(1,1)+1,1)=x(isearch+search_area*(isearch_time-1))+xmin-1;
                    Sreplace2(Sreplace2(1,1)+1,2)=y(isearch+search_area*(isearch_time-1))+ymin-1;
                    Sreplace2(Sreplace2(1,1)+1,3)=z(isearch+search_area*(isearch_time-1))+zmin-1;
                end
                
                
                if matching_number>1
                    Sadd2{matching_number-1}(1,1)=Sadd2{matching_number-1}(1,1)+1;
                    Sadd2{matching_number-1}(Sadd2{matching_number-1}(1,1)+1,1)=x(isearch+search_area*(isearch_time-1))+xmin-1;
                    Sadd2{matching_number-1}(Sadd2{matching_number-1}(1,1)+1,2)=y(isearch+search_area*(isearch_time-1))+ymin-1;
                    Sadd2{matching_number-1}(Sadd2{matching_number-1}(1,1)+1,3)=z(isearch+search_area*(isearch_time-1))+zmin-1;
                end
                
            end
        end
    end
    % Now all the matching is done, calculation the rate of enlargement and record the number
    Sreplace2(1,3)=Sreplace2(1,1)/Sreplace(1,1);
    if nofnewcell>0
        for ilabel=1:nofnewcell
            Sadd2{ilabel}(1,3)=Sadd2{ilabel}(1,1)/Sadd{ilabel}(1,1);
        end
    end
    Sadd=Sadd2;Sreplace=Sreplace2;
else
    % If nofnewcell == 0, the cluster is a single cell
    SS(1,3)=SS(1,1)/Sreplace(1,1);
    Sreplace=SS;
end




