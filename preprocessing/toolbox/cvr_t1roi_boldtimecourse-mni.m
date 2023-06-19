function [sum1, sum2, sum3,name1,name2,name3,nvox1,nvox2,nvox3] = cvr_t1roi_boldtimecourse(path_mpr,name_mpr,boldvols,rmskfile,outpath)
% calculate the CBF ROI average for three segmentation levels (Type1-L5, Type1-L3, and Type1-L2)
% yli20160627

% read the mutilevel_lookup_table.txt to determine the version of T1
% multiatlas and ROI lookups for all levels
roi_lookup_file = [path_mpr filesep 'multilevel_lookup_table.txt'];
roi_lookup_all  = read_roi_lookup_tabl(roi_lookup_file);

label_idx = str2num(char(roi_lookup_all{1,1}));
label_num = length(label_idx);
atlas_ver = [num2str(label_num)];


% read the roi name lists for the three segmentation levels
roitypes = {'Type1-L2';
            'Type1-L3';
            'Type1-L5'}; % order: type1-L2, type1-L3, type1-L5
roi_lookup_tbl = {roi_lookup_all{1,5};
                  roi_lookup_all{1,4};
                  roi_lookup_all{1,2}}; % orduer: type1-L2, type1-L3, type1-L5
              
roi_stats_file = spm_select('FPList',path_mpr,['^' name_mpr '.*' atlas_ver '.*MNI_stats.txt$']);
roi_stats_file = roi_stats_file(1,:);
roi_lists_info = read_roi_lists_info(roi_stats_file,roitypes);
roi_lists_info{3,2} = label_num;
roi_lists_info{3,3} = roi_lookup_all{1,2};


% get mask and CBF maps
P0 = spm_select('FPList',path_mpr,['^' name_mpr '.*' atlas_ver '.*MNI_2mm.img$']); % apply 2mm mni labeled mask
roimaskfile = P0;

V0      = spm_vol(roimaskfile);
V1      = spm_vol(boldvols);
V3      = spm_vol(rmskfile);
mskv    = spm_read_vols(V0);
bold    = spm_read_vols(V1); bold(isnan(bold)) = 0;
msk1    = spm_read_vols(V3); msk1(isnan(msk1)) = 0;

matsize = size(bold);
boldvec = reshape(bold, [],matsize(4));

sum1 = zeros(length(roi_lists_info{1,3}),matsize(4));
sum2 = zeros(length(roi_lists_info{2,3}),matsize(4));
sum3 = zeros(length(roi_lists_info{3,3}),matsize(4));
name1 = cell(length(roi_lists_info{1,3}),1);
name2 = cell(length(roi_lists_info{2,3}),1);
name3 = cell(length(roi_lists_info{3,3}),1);
nvox1 = zeros(length(roi_lists_info{1,3}),1);
nvox2 = zeros(length(roi_lists_info{2,3}),1);
nvox3 = zeros(length(roi_lists_info{3,3}),1);

% Combine parcellations to Type1-L2 or Type1-L3 segmentations
for tt = 1:2
    roi_tbl = roi_lookup_tbl{tt};
    roi_lst = roi_lists_info{tt,3};
    
    % label the mask with other level of segmentation
    seg_num = length(roi_lst);
    seg_idx = zeros(size(roi_lst));
    segmask = zeros(size(mskv));
    
    for ii = 1:seg_num
        seg_idx(ii) = ii;
        for kk = 1:label_num
            if strcmp(roi_lst{ii},roi_tbl{kk})
                segmask(mskv == label_idx(kk)) = ii;
            end
        end
    end
    
%     % write out the segmented mask volume
%     outVol = V0;
%     outVol.fname = [outpath filesep name_mpr '_' num2str(seg_num) 'segments.img'];
%     spm_write_vol(outVol, segmask);
    
    % calculate the ROI CBF average and print to text file
    for ii = 1:seg_num
        idxvox_seg = logical( (segmask == seg_idx(ii)) .* (msk1 > 0.5) ); % count vox only within ASL mask (in case of ASL partial coverage)
        seg_bold = mean( boldvec( idxvox_seg(:),:), 1);
        seg_nvox = length( find( idxvox_seg == 1 ) );
        seg_name = roi_lst{ii};
        if tt == 1
            sum1(ii,:) = seg_bold;
            name1{ii}  = seg_name;
            nvox1(ii)  = seg_nvox;
        elseif tt == 2
            sum2(ii,:) = seg_bold;
            name2{ii}  = seg_name;
            nvox2(ii)  = seg_nvox;
        end
    end
end

% % calculate the ROI CBF average of Type1-L5 segmentations and print
% for ii = 1:label_num
%     idxvox_seg = logical( (mskv == label_idx(ii)) .* (msk1 > 0.5) );
%     seg_bold = mean( boldvec( idxvox_seg(:),:), 1);
%     seg_nvox = length(bold(idxvox_seg));
%     seg_name = roi_lists_info{3,3}{ii};
%     sum3(ii,:) = seg_bold;
%     name3{ii}  = seg_name;
%     nvox3(ii)  = seg_nvox;
% end

end


function roi_lookup_tabl = read_roi_lookup_tabl(roi_lookup_file)
fileid = fopen(roi_lookup_file);
title  = textscan(fileid,'%s',1,...
                    'delimiter','\n');

roi_lookup_tabl = textscan(fileid,repmat('%s ',1,11),...
                    'delimiter',{' ','\b','\t'},'MultipleDelimsAsOne',1);
                
fclose(fileid);
end


function roi_lists_info = read_roi_lists_info(roi_stats_file,roitypes)
% return the level tag, # of rois, & roi list for each segmentation level

fileid  = fopen(roi_stats_file);
tmp     = textscan(fileid,'%s','delimiter','\n','whitespace','');
fclose(fileid);

alllines = tmp{:};
nline    = length(alllines);
ntype    = length(roitypes);

tmpinfo  = cell(ntype,3);
iline = 1;
while iline < nline
    splitStr = regexp(alllines{iline,1},'[ \t\b]','split');
    
    for itype = 1:ntype
        if strcmp(splitStr{1},roitypes{itype})
            tmpinfo{itype,1} = splitStr{1};
            
            roilist     = cell(1,1);
            roicount    = 0;
            iline       = iline + 2;
            while ~isempty(regexp(alllines{iline,1},'.img','ONCE'))
                roicount = roicount+1;
                
                splitStr2 = regexp(alllines{iline,1},'[ \t]','split');
                roilist{roicount,1} = splitStr2{2};
                
                iline = iline + 1;
            end
            
            tmpinfo{itype,2} = roicount;
            tmpinfo{itype,3} = roilist;
        end
    end
    
    iline = iline + 1;
end

roi_lists_info = tmpinfo;
end

