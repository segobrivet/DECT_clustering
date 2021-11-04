
%% LOAD DATASET AND GROUND TRUTH SEGMENTATION


%%%%%%%%%%%%%%%

close all;
clear; clc;

addpath('utils');
addpath('utils_organs');

data_folder = "C:\\Users\\Segolene\\Documents\\Canada\\McGill\\PhD\\Multi-energy CT\\data\\SRG_MultiEnergy";
segm_folder = "C:\\Users\\Segolene\\Documents\\Canada\\McGill\\PhD\\Multi-energy CT\\data\\SRG_MultiEnergy";
% data_folder = "/Users/Shared/datasts/HNSCC/Multi-energy";
% segm_folder = "/Users/Shared/datasts/HNSCC/Multi-energy";

% patient_names = ["SRG12_MultiEnergy","SRG14_MultiEnergy","SRG15_MultiEnergy","SRG18_MultiEnergy","SRG21_MultiEnergy","SRG23_MultiEnergy","SRG24_MultiEnergy","SRG32_MultiEnergy","SRG92_MultiEnergy","SRG94_MultiEnergy"];
patient_names = ["SRG3_MultiEnergy","SRG7_MultiEnergy","SRG8_MultiEnergy","SRG14_MultiEnergy"];
segm_type = 'organs';
additional_vars.organ_id = {{3},{4}};
additional_vars.verbose = 1;


for patient_name = patient_names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



[char_subj, subject_nb, subject_id, org_str_list] = get_subject_names(patient_name, segm_type, additional_vars.organ_id);
study_type = '';
medic_img_struct = '';

%%% Build Ground truth segmentation
gt_path_name_list = {};
build_subject = 0;
for s=1:length(org_str_list)
    gt_path_name_list{s} = fullfile('..','data_organs',[subject_id,'_',org_str_list{s},'_GT.mat']);
    if ~isfile(gt_path_name_list{s}), build_subject = 1; end
end

segm_vol_full_list = {};
if build_subject == 0
    % load all segmentations from saved variables
    for c=1:length(gt_path_name_list)
        segm_vol_full_list{c} = load(gt_path_name_list{c}).segm_vol_full;
        if additional_vars.verbose
            disp("Ground truth '" + gt_path_name_list{c} + "' has been loaded."); disp(' ');
        end
    end
    
else
    % import and build segmentations
    [medic_img_struct, ~, segm_vol_full_list, study_type] = ...
        load_organ_segm(segm_folder, patient_name, segm_type, additional_vars.organ_id);

    % save segmentations
    for c=1:length(segm_vol_full_list)
        segm_vol_full = segm_vol_full_list{c};
        org_str = org_str_list{c};
        save(gt_path_name_list{c},'segm_vol_full');
        if additional_vars.verbose
            disp("Ground truth saved with the name '" + gt_path_name_list{c} + "'."); disp(' ');
        end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%% Build subject scan sections
subject_path_name_list = {};
build_subject = 0;
for s=1:length(org_str_list)
    subject_path_name_list{s} = fullfile('..','data_organs',[subject_id,'_',org_str_list{s},'.mat']);
    if ~isfile(subject_path_name_list{s}), build_subject = 1; end
end

subject_list = {};
slices_range_list = {};
if build_subject == 0
    % load all subject scan sections from saved variables
    for c=1:length(subject_path_name_list)
        subject_list{c} = load(subject_path_name_list{c}).subject;
        if additional_vars.verbose
            disp("Scan section '" + subject_path_name_list{c} + "' has been loaded."); disp(' ');
        end
        
        % find slices to work on
        lin_ind = find(segm_vol_full_list{c});
        [~,~,z] = ind2sub(size(segm_vol_full_list{c}),lin_ind);
        slices_range_list{c} = (min(z):max(z))'; % in Matlab indices (then subtract 1 for Slicer indices)
    end
    
else
    % import and build subject full scan
    
%     % in case GT was saved but not subject... the study type would be unkown
%     if ~exist('study_type','var')     % 'study_type' var exists already in this script
%         [~, ~, ~, study_type] = load_organ_segm(segm_folder, patient_name, segm_type, additional_vars.organ_id);
%     end

    subject_full = build_patient_img_srg(data_folder, patient_name, study_type);
    
    % if GT and subject don't have the same nb of slices,
    % either duplicate some GT slices to fit subject slices, or select some slices in GT to match the fewer in subject
    if size(segm_vol_full_list{1},3) ~= size(subject_full{1},3)
        % select_ind = adjust_GT_slices(size(segm_vol_full_list{1},3), size(subject_full{1},3));
        [select_ind,vrb] = adjust_GT_subj_slices(size(segm_vol_full_list{1},3), size(subject_full{1},3));
        
        switch(vrb)
            case 'subj'
                for kev=1:length(subject_full)
                    subject_full{kev} = subject_full{kev}(:,:,select_ind);
                end
                    
            case 'gt'      
                for c=1:length(segm_vol_full_list)
                    segm_vol_full = segm_vol_full_list{c}(:,:,select_ind);
                    segm_vol_full_list{c} = segm_vol_full;
                    
                    % rewrite ground truth
                    save(gt_path_name_list{c},'segm_vol_full');
                    if additional_vars.verbose
                        disp("Ground truth '" + gt_path_name_list{c} + "' has been updated."); disp(' ');
                    end
                end
            otherwise
                % never happens
        end
    end

    for c=1:length(subject_path_name_list)
        % find slices to work on
        lin_ind = find(segm_vol_full_list{c});
        [~,~,z] = ind2sub(size(segm_vol_full_list{c}),lin_ind);
        slices_range_list{c} = (min(z):max(z))'; % in Matlab indices (then subtract 1 for Slicer indices)

        subject = cell(1, length(subject_full));
        for kev=1:length(subject_full)
            subject{kev} = subject_full{kev}(:,:,slices_range_list{c});   % slices_range-1 ?
        end
%         subject_list{c} = subject;
        % save subject cut
        save(subject_path_name_list{c},'subject');
        if additional_vars.verbose
            disp("Scan section saved with the name '" + subject_path_name_list{c} + "'."); disp(' ');
        end
    end
end


end