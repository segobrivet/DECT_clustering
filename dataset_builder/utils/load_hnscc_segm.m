
function [segm_vol_full, medic_img_struct]  = load_hnscc_segm(segm_folder, patient_name, segm_type)
% find the desired segmentation file based on folder names

%%% ARGUMENTS
% data_folder           : string containing the data folder path
% patient_name          : string containing the patient identification name (ex: "HNSCC5")
% segm_type             : string containing the segmentation type

% find segmentation file
switch segm_type
    case 'tumor'
        pn = char(patient_name);
        tumor_folder = fullfile(segm_folder,patient_name,'ROI mask', ['Hnscc',pn(6:end)]);
    otherwise
        error("Only tumor segmentations exist for HNSCC.");
end

% segmentation volume
segm_vol_full = read_DICOM_from_dir(tumor_folder, []);
segm_vol_full = permute(segm_vol_full, [2 1 3]);

% segmentation header
dname=dir(fullfile(tumor_folder,'IM*'));
medic_img_struct = dicominfo(fullfile(tumor_folder,dname(1).name));



% % % % % % % % % %%%%%  REFINEMENTS IF TOO BIG  %%%%% % % % % % % % % % %

%%% If tumor is in several parts with discontinuities on Z dim,
%%% select only the biggest part of Z slices
[~,~,slic_z] = ind2sub(size(segm_vol_full),find(segm_vol_full));
slic_z = unique(slic_z);
ds = diff(slic_z);
ch = find(ds-1);  % indices of discontinuities
if ~isempty(ch)
    ch = [0; ch; length(slic_z)];
    [~,ind] = max(ch(2:end)-ch(1:end-1));
    selected_slices = slic_z((ch(ind)+1):ch(ind+1));
    slic_z = selected_slices;
end


%%%% if tumor is too big, take only middle slices

if length(slic_z) > 20
    mid = ceil(length(slic_z)/2);
    selected_slices = slic_z(mid-10 : mid+9); % 20 middles slices
end


%%% Replace by 0 all other parts of tumor if necessary

if exist('selected_slices','var')
    
    segm_vol = segm_vol_full(:,:,selected_slices);
    segm_vol_full = zeros(size(segm_vol_full)); 
    segm_vol_full(:,:,selected_slices) = segm_vol;
end

% TODO: to check

end
