% Load data and build dataset for ground truth and curves and prepare plots

function res = build_patient_data(machine_type, patient_nm, take_all_slices, org_ids, roi_radius, max_slic_subplot, show_slices)

    switch(machine_type)
        case('sego_laptop')
            if strcmp(patient_nm(1:4),'subj') && (~strcmp(patient_nm(end-4:end),'tumor'))
                load(['../../data_organs/',patient_nm,'_organs_',org_ids,'_GT.mat']);
                load(['../../data_organs/',patient_nm,'_organs_',org_ids,'.mat']);
            else
                load(['../../data/',patient_nm,'_GT.mat']);
                load(['../../data/',patient_nm,'.mat']);
            end
        case('sego_clean')
            if strcmp(patient_nm(1:4),'subj') && (~strcmp(patient_nm(end-4:end),'tumor'))
                load(['../data_organs/',patient_nm,'_organs_',org_ids,'_GT.mat']);
                load(['../data_organs/',patient_nm,'_organs_',org_ids,'.mat']);
            else
                load(['../data_tumor/',patient_nm,'_GT.mat']);
                load(['../data_tumor/',patient_nm,'.mat']);
            end
        case('Peter')
            load(['../data/',patient_nm,'_GT.mat']);
            load(['../data/',patient_nm,'.mat']);
        case('GPU')
            if strcmp(patient_nm(1:4),'subj') && (~strcmp(patient_nm(end-4:end),'tumor'))
                load(['../data/matfiles/',patient_nm,'_organs_',org_ids,'_GT.mat']);
                load(['../data/matfiles/',patient_nm,'_organs_',org_ids,'.mat']);
            else
                load(['../data/matfiles/',patient_nm,'_GT.mat']);
                load(['../data/matfiles/',patient_nm,'.mat']);
            end
		case('GPU_clean')
            if strcmp(patient_nm(1:4),'subj') && (~strcmp(patient_nm(end-4:end),'tumor'))
                % ?
            else
                load(['../../data/data_tumor/',patient_nm,'_GT.mat']);
                load(['../../data/data_tumor/',patient_nm,'.mat']);
            end
        case('Faicel')
            ...
        otherwise
            error('unknown machine type');
    end
    
    % Get tumor slices
    [~,~,slic_min] = ind2sub(size(segm_vol_full),find(segm_vol_full,1,'first')); % lower slice containing a tumor
    [~,~,slic_max] = ind2sub(size(segm_vol_full),find(segm_vol_full,1,'last')); % higher slice containing a tumor
    % Selecting all tumor slices
    if take_all_slices   % take all tumor slices
        slic_inds = 1:(slic_max-slic_min+1);
    else   % take 8 middle slices
        slic_inds = round((slic_max-slic_min+1)/2)-4 : round((slic_max-slic_min+1)/2)+3;
        if min(slic_inds) <= 0, slic_inds = slic_inds-slic_inds(1)+1; end
        if max(slic_inds) > (slic_max-slic_min+1), slic_inds = slic_inds(slic_inds<=(slic_max-slic_min+1)); end
    end
    
    % Selecting the 6 slices in the middle of the tumor to plot on figure
%     max_slic_subplot = 6;
    if length(slic_inds) >= max_slic_subplot
        switch(show_slices)
            case 'middle'
                slic_show = round(length(slic_inds)/2)-floor((max_slic_subplot-1)/2) : round(length(slic_inds)/2)+floor(max_slic_subplot/2);
            case 'top'
                slic_show = (length(slic_inds)-(max_slic_subplot-1)):length(slic_inds);
            case 'bottom'
                slic_show = 1:max_slic_subplot;
        end
    else
        slic_show = 1:length(slic_inds);
    end
    
    
    % contour of body
    [subj_slic, focus, xy_min, xy_max] = mask_outOfBody(permute(subject{1}(:,:,slic_inds),[2 1 3]));
    
    
    % Select a ROI of the 3D volume    
    gr_truth = segm_vol_full(:,:,slic_min+slic_inds-1);
    s = regionprops(gr_truth,'centroid');
    try
        assert(length(s.Centroid) == 3);
    catch 
        error("Ground truth tumor contains more than 1 connected component");
    end
    tc = round(s.Centroid(1)); tr = round(s.Centroid(2));
    
    rmax = tr + floor(roi_radius/2);
    rmin = tr - (ceil(roi_radius/2)-1);
    if rmax > xy_max(1)
        shift = rmax - xy_max(1);
        rmax = xy_max(1);
        rmin = rmin - shift;
        if rmin < xy_min(1)
            rmin = xy_min(1);
        end
    end
    cmax = tc + floor(roi_radius/2);
    cmin = tc - (ceil(roi_radius/2)-1);
    if cmax > xy_max(2)
        shift = cmax - xy_max(2);
        cmax = xy_max(2);
        cmin = cmin - shift;
        if cmin < xy_min(2)
            cmin = xy_min(2);
        end
    end
    
    
    % Remove air voxel
    gr_truth_dil = zeros(size(gr_truth));
    for i=1:length(slic_inds)
        gr_truth_dil(rmin:rmax,cmin:cmax,i) = ones(rmax-rmin+1,cmax-cmin+1);  % for a square
%         gr_truth_dil(:,:,i) = imdilate(gr_truth(:,:,i),strel('disk',roi_radius)); % for a circle
    end
    lin_obj = find(gr_truth_dil);  % ROI around tumor;
    
    % Remove air voxel
    lin_air = subj_slic(lin_obj)<-500;
    lin_obj = lin_obj(~lin_air);
    
    % find coordinates of ROIs and frame of the plot
    [row_obj, col_obj, z_obj] = ind2sub(size(gr_truth),lin_obj);
    coord = [row_obj, col_obj, z_obj+slic_inds(1)-1];
%     % following 2 lines for a circle
%     rmin = max(min(row_obj)-10,1); cmin = max(min(col_obj)-10,1);
%     rmax = min(max(row_obj)+10,size(subj_slic,1)); cmax = min(max(col_obj)+10,size(subj_slic,2));
    
    % Build tumor class label map
    lin_tum = find(gr_truth);  % tumor;
    klas_tum = ismember(lin_obj,lin_tum);
    
    %% Select curves of voxels in ROI
    
    % line 'i' stores the corresponding coordinates and curves
    
    decay_curves = zeros(size(coord,1),21);
    for r=1:size(coord,1)
        for kev=1:21
            decay_curves(r,kev) = subject{kev}(coord(r,2),coord(r,1),coord(r,3));  % careful to transposed subject: dim 2 <-> 1.
        end
    end
    Y = decay_curves;
    Y = zscore(Y);%Y - ones(length(Y),1)*mean(Y,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    res.Y = Y;
    res.decay_curves = decay_curves;
    res.klas_tum = klas_tum;
    res.lin_obj = lin_obj;
    res.gr_truth = gr_truth;
    res.coord = coord;
    res.rmin = rmin;
    res.rmax = rmax;
    res.cmin = cmin;
    res.cmax = cmax;
    res.slic_show = slic_show;
    res.slic_inds = slic_inds;
    res.subj_slic = subj_slic;
    res.slic_min = slic_min;
    res.slic_min_idx = slic_inds(1)+slic_min-2;
    res.focus = focus;
    
    end
    