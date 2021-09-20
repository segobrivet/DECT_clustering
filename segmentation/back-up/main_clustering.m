%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Learning spatial mixture of functional regression models for spectral image clustering
%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Choose specifications

clear all; close all; 
clc;

addpath('utils');

%%%%%%%%%%%%%%%%%%%%%%%%%%% MACHINE and FOLDER NAMES %%%%%%%%%%%%%%%%%%%%%%
% machine_type = 'sego_laptop';
machine_type = 'sego_clean';
% machine_type = 'generic_data';
% machine_type = 'Peter';
% machine_type = 'GPU';

% results_folder_name = 'results_generic_data';
% results_folder_name = 'results_clustering_L07_K60';
% results_folder_name = 'results_small';
results_folder_name = 'results_test';
% results_folder_name = 'results_organs_Fun_R08_SRG21/fixSpatKmeans';
% results_folder_name = 'results_3d_bigROIs_Fun_R07';
% results_folder_name = 'results_3d_bigROIs_Fun_R08_iCurves_gSpat';
% results_folder_name = 'results_3d_bigROIs_extraFeat_R09_cl70';
% results_folder_name = 'results_stats';


%%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
strategy = 'FunSRM'; init_kmeans_dep = 1;
% strategy = 'ExtraFeat-GMM';
% strategy = 'kmeans';
% strategy = 'imsegkmeans';

% mixingOption = 'softmax';
mixingOption = 'gaussian';

lambda = 0.07;

% model = "PRM"; % Polynomial regression mixture
% model = "SRM"; % Spline regression mixture
model = "bSRM";% B-Spline regression mixture

% p = 3; % polynomial regression degree
% spline_order = 4; % 2: linear, 3: qudratic, 4: cubic, etc, for (b-)spline spatial regression mixture
Bspline_order = 4; % 2: linear, 3: qudratic, 4: cubic, etc, for (b-)spline spatial regression mixture
nknots = 10; % fixed number of internal knots, for (b-)spline spatial regression mixture


%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROI size
roi_radius = 130;

% tissue enhancement window
lvl = 150;  % 40  % 70
wdw = 700;  % 350  % 500

org_ids = '1-2';

%%%%%%%%%%%%%%%%%%%%%%%%%%% CHOOSE PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
show_slices = 'middle';
% show_slices = 'bottom';
% show_slices = 'top';
max_slic_subplot = 6;

plot_clusterCurves = 0;

plot_initResults = 1;
plot_filter3 = 0;
plot_filter5 = 0;
plot_initMerge = 0;
plot_mergeFilter3 = 1;
plot_rab = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%% TO LOAD A SAVED MODEL %%%%%%%%%%%%%%%%%%%%%%%%%
load_saved_mdl = 1;
pat = [];
dic = [];
jac = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%% NB OF RESTART %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nbr_EM_runs = 1;

%% end of specifications

for K = 60   % 300 %[40,50,60] % number of clusters

%% To read saved result variables in a folder and plot the images
% fls = dir(fullfile(results_folder_name,'/subj*stats_red.mat'));
% for nm = 1:length(fls)
%     mixstatspath = fullfile(results_folder_name,fls(nm).name);
%     str = split(fls(nm).name,'_');
%     patient_name = 'subject8_tumor'; rdn = str{3};
%     roi_radius = str{5}; roi_radius = str2num(roi_radius(4:end));
%     K = str{6}; K = str2num(K(3:end));

    % % OR % %
    
% % fls = dir(fullfile(results_folder_name,'/subject21_127*mixstats_red.mat'));
% fls = dir(fullfile(results_folder_name,'/HNSCC18*cl50_mixstats_red.mat'));
% for nm = 1:length(fls)
% %     imsegkmeanspath = fullfile(results_folder_name,fls(nm).name);
%     mixstatspath = fullfile(results_folder_name,fls(nm).name);
% %     mixmodelpath = [mixstatspath(1:end-12),'mixmodel.mat'];  % (1:end-17)
%     str = split(fls(nm).name,'_');
%     patient_name = str{1}; rdn = str{2};
%     roi_radius = str{4}; roi_radius = str2num(roi_radius(4:end));
%     K = str{5}; K = str2num(K(3:end));


% % for visualising init MyKmeans
% % fls = dir(fullfile(results_folder_name,'/solInitKmeans*.mat'));
% fls = dir('solInitKmeansSpat*.mat');
% for nm = 1:length(fls)
% %     mixstatspath = fullfile(results_folder_name,fls(nm).name);
%     mixstatspath = fullfile(fls(nm).name);
%     str = split(fls(nm).name,'_');
%     patient_name = 'subject21';
%     roi_radius = 90;
%     K = str{2}; K = str2num(K(3:end));
%     rdn = str{3}; rdn = rdn(1:end-4);

%%

% for patient_name = ["subject3","subject12","subject14","subject18","subject21","subject23","subject7"]

for patient_name = ["HNSCC1"]
    
% for patient_name = ["subject8_tumor","HNSCC2","HNSCC3","HNSCC5","HNSCC8","HNSCC9","HNSCC10","HNSCC11","HNSCC12","HNSCC13","HNSCC15","HNSCC17","HNSCC18","HNSCC26"]
    
    
    %% Data (spectral image)

    close all; 
    
    patient_nm = char(patient_name); disp(patient_nm);
    % fn_save_pdf = '';  % if don't want to save fig as pdf
    fn_save_pdf = fullfile(results_folder_name, patient_nm);  % to save fig as pdf
    
    %%% Build dataset for a specific patient
    if strcmp(machine_type,'generic_data')
        load('../generic_data/subject.mat')
        load('../generic_data/gr_truth.mat')
        
        slic_inds = 1:size(gr_truth,3);
        gr_truth_dil = zeros(size(gr_truth));
        for i=1:size(gr_truth,3)
            gr_truth_dil(:,:,i) = imdilate(gr_truth(:,:,i),strel('disk',roi_radius));
        end
        lin_obj = find(gr_truth_dil);  % ROI around tumor;
        % Remove air voxel
        [subj_slic, focus, xy_min, xy_max] = mask_outOfBody(subject{1}(:,:,slic_inds));
        lin_air = subj_slic(lin_obj)<-500;
        lin_obj = lin_obj(~lin_air);
        % find coordinates of ROIs and frame of the plot
        [row_obj, col_obj, z_obj] = ind2sub(size(gr_truth),lin_obj);
        coord = [row_obj, col_obj, z_obj];
        rmin = max(min(row_obj)-10,1); cmin = max(min(col_obj)-10,1);
        rmax = min(max(row_obj)+10,size(subj_slic,1)); cmax = min(max(col_obj)+10,size(subj_slic,2));
        
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
%         Y = zscore(Y);%Y - ones(length(Y),1)*mean(Y,1);
    
        slic_show = 1:6;
        slic_min_idx = 1;
        slic_min = 1;
        
    else
        take_all_slices = 0;
        res = build_patient_data(machine_type, patient_nm, take_all_slices, org_ids, roi_radius, max_slic_subplot, show_slices);

        Y = res.Y; decay_curves = res.decay_curves; gr_truth = res.gr_truth;
        klas_tum = res.klas_tum; lin_obj = res.lin_obj; subj_slic = res.subj_slic;
        slic_show = res.slic_show; slic_min_idx = res.slic_min_idx; slic_min = res.slic_min;
        rmin = res.rmin; rmax = res.rmax; cmin = res.cmin; cmax = res.cmax; coord = res.coord;
    end
    
    [n, m] = size(Y);
    V(:,3) = coord(:,3)*2;   % spacing in Z is 1.25mm, spacing in X and in Y is 0.61mm, so Z is 2*bigger than X-Y.
    T = linspace(0, 1, 21);  % T = linspace(40,140,21); %1:21;
    
    Curves.spatialcoord = V;
    Curves.abscissas = T;
    Curves.ordinates =  Y;
    

    %% Plot original image with tissue enhancement
    
    subj_slic(subj_slic<(lvl-(wdw/2))) = lvl-(wdw/2);
    subj_slic(subj_slic>(lvl+(wdw/2))) = lvl+(wdw/2);
    
    len_sp = min(ceil(length(slic_show)/2),ceil(max_slic_subplot/2));
    slics = subj_slic(rmin:rmax,cmin:cmax,:);
    
    for i=1:length(slic_show)
%         tumor_contour_list{i} = bwboundaries(segm_vol_full(:,:,slic_inds(1)+slic_show(i)-1+slic_min-1 -1));   % shift of +1 in segm_vol_full/subject?
        tumor_contour_list{i} = bwboundaries(gr_truth(:,:,slic_show(i)));
    end
    
    % original clustering
    if plot_initResults
        fig_slic = figure('units','normalized','outerposition',[0 0 1 1]);
        for i=1:length(slic_show)
            if i > max_slic_subplot, break; end
            subplot(len_sp,4,2*i-1);
            imshow(slics(:,:,slic_show(i)),[]), hold on, title("slice " + num2str(slic_inds(1)+slic_show(i)-1+slic_min-2)) % index showed on 3DSlicer
            plot_tumor_contour(tumor_contour_list{i}, [rmin, cmin], [0,0.5,1]);
        end
    end
    
    % original clustering with filter
    if plot_filter3
        fig_slic2 = figure('units','normalized','outerposition',[0 0 1 1]); 
        for i=1:length(slic_show)
            if i > max_slic_subplot, break; end
            subplot(len_sp,4,2*i-1);
            imshow(slics(:,:,slic_show(i)),[]), hold on, title("slice " + num2str(slic_inds(1)+slic_show(i)-1+slic_min-2)) % index showed on 3DSlicer
            plot_tumor_contour(tumor_contour_list{i}, [rmin, cmin], [0,0.5,1]);
        end
    end
    
    % original clustering with other filter
    if plot_filter5
        fig_slic3 = figure('units','normalized','outerposition',[0 0 1 1]); 
        for i=1:length(slic_show)
            if i > max_slic_subplot, break; end
            subplot(len_sp,4,2*i-1);
            imshow(slics(:,:,slic_show(i)),[]), hold on, title("slice " + num2str(slic_inds(1)+slic_show(i)-1+slic_min-2)) % index showed on 3DSlicer
            plot_tumor_contour(tumor_contour_list{i}, [rmin, cmin], [0,0.5,1]);
        end
    end
    
    % original clustering with merging step
    if plot_initMerge
        fig_slic4 = figure('units','normalized','outerposition',[0 0 1 1]);  
        for i=1:length(slic_show)
            if i > max_slic_subplot, break; end
            subplot(len_sp,4,2*i-1);
            imshow(slics(:,:,slic_show(i)),[]), hold on, title("slice " + num2str(slic_inds(1)+slic_show(i)-1+slic_min-2)) % index showed on 3DSlicer
            plot_tumor_contour(tumor_contour_list{i}, [rmin, cmin], [0,0.5,1]);
        end
    end
    
    % original clustering with merging step and filter
    if plot_mergeFilter3
        fig_slic5 = figure('units','normalized','outerposition',[0 0 1 1]);  
        for i=1:length(slic_show)
            if i > max_slic_subplot, break; end
            subplot(len_sp,4,2*i-1);
            imshow(slics(:,:,slic_show(i)),[]), hold on, title("slice " + num2str(slic_inds(1)+slic_show(i)-1+slic_min-2)) % index showed on 3DSlicer
            plot_tumor_contour(tumor_contour_list{i}, [rmin, cmin], [0,0.5,1]);
        end
    end
    
    if plot_rab   % if need another plot
        fig_slic6 = figure('units','normalized','outerposition',[0 0 1 1]);  
        for i=1:length(slic_show)
            if i > max_slic_subplot, break; end
            subplot(len_sp,4,2*i-1);
            imshow(slics(:,:,slic_show(i)),[]), hold on, title("slice " + num2str(slic_inds(1)+slic_show(i)-1+slic_min-2)) % index showed on 3DSlicer
            plot_tumor_contour(tumor_contour_list{i}, [rmin, cmin], [0,0.5,1]);
        end
    end
    
    
    
    %% model specification
    switch(model)
        case('PRM')
            regressionOptions.basis = 'polynomial';
            regressionOptions.p = p;
        case 'SRM'
            regressionOptions.basis = 'spline';
            regressionOptions.spline_order = spline_order;
            regressionOptions.nknots = nknots;
        case('bSRM')
            regressionOptions.basis = 'B-spline';
            regressionOptions.Bspline_order = Bspline_order;
            regressionOptions.nknots = nknots;
        otherwise
            error('unknown model type');
    end
    %%
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            MODEL FITTING                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% SRM MODEL FITTING

    if load_saved_mdl

        switch(strategy)
            case('imsegkmeans')
                img_4D = zeros(rmax-rmin+1, cmax-cmin+1, length(slic_inds), 3);
                kev_idx = 1;
                for kev=[1,11,21]
                    img_4D(:,:,:,kev_idx) = single(permute(map_to_0_1(subject{kev}(cmin:cmax,rmin:rmax,slic_inds)), [2 1 3]));  % .*focus
                    kev_idx = kev_idx+1;
                end
                
                load(imsegkmeanspath);
                
            case('kmeans')
                load(mixstatspath);
                mixstats.klas = sol_km.klas;
                mixstats.Muk = sol_km.muk;
                
            otherwise
%                 load(mixmodelpath);
%                 load(mixstatspath);

                load(mixstatspath);
                mixstats.klas = mixstats_red.klas;
                mixstats.Muk = mixstats_red.Muk;
        end

    else

        disp('Starting EM...')
        tic;
        switch(strategy)

            case('FunSRM')
                
                nb_EM = 1; init_kmeans = init_kmeans_dep;
                while nb_EM > 0 && nb_EM < 7
                    try
                        [mixmodel, mixstats] = learn_SRM_EM(Curves, K, mixingOption, regressionOptions, nbr_EM_runs, lambda, init_kmeans);
                        disp(['Optim successful after ',num2str(nb_EM),' EM run(s).']);
                        nb_EM = 0;
                    catch
                        disp("Optim failed, try again.");
                        nb_EM = nb_EM + 1;
                        init_kmeans = 0;
                    end
                end
                if nb_EM == 7
                    disp("Switch to next patient.");
                    continue;
                end

            case('ExtraFeat-GMM')
                % construct B splines features
                [~, B] = designSRM(Curves.spatialcoord, Curves.abscissas, mixingOption, regressionOptions);
                dimBeta = size(B,2);
                Betas = zeros(n, dimBeta);
%                     figure;
                for i=1:n
                    betai = (B'*B)\B'*Y(i,:)';
                    Betas(i,:) = betai;
%                         plot(Curves.abscissas, Curves.ordinates(i,:),'o');
%                         hold on, plot(Curves.abscissas, B*betai,'r');
%                         pause
%                         clf
                end
                %
                Betas = zscore(Betas);
                Curves.coefficients = Betas;

                nb_EM = 1;
                while nb_EM > 0 && nb_EM < 7
                    try
                        [mixmodel, mixstats] = learn_SRM_EM_Gauss(Curves, K, mixingOption, regressionOptions, nbr_EM_runs, lambda);
                        disp(['Optim successful after ',num2str(nb_EM),' EM run(s).']);
                        nb_EM = 0;
                    catch
                        disp("Optim failed, try again.");
                        nb_EM = nb_EM + 1;
                    end
                end
                if nb_EM == 7
                    disp("Switch to next patient.");
                    break;
                end
                
                mixstats.Muk = mixmodel.Muk';

            case('kmeans')
                [mixstats.klas, Muk] = kmeans(Curves.ordinates, K, 'MaxIter', 500, 'Display', 'iter');
                mixstats.Muk = Muk';
%                 max_iter_kmeans = 500;
%                 n_tries_kmeans = 1; %3
%                 verbose_kmeans = 1;
%                 sol_km = myKmeans(Curves.ordinates, K , n_tries_kmeans, max_iter_kmeans, verbose_kmeans);
%                 mixstats.klas = sol_km.klas;
%                 mixstats.Muk = sol_km.muk;
% %                 mixmodel.Muk = sol_km.muk;

            case('imsegkmeans')

                img_4D = zeros(rmax-rmin+1, cmax-cmin+1, length(slic_inds), 3);
                kev_idx = 1;
                for kev=[1,11,21]
                    img_4D(:,:,:,kev_idx) = single(permute(map_to_0_1(subject{kev}(cmin:cmax,rmin:rmax,slic_inds)), [2 1 3]));  % .*focus
                    kev_idx = kev_idx+1;
                end
                [X,Y,Z] = meshgrid(1:size(img_4D,2),1:size(img_4D,1),(1:size(img_4D,3)));  % Z *1 or *5 (coord spacing) gives equal results
                % featureSet = cat(4,img_4D, X, Y, Z);
                featureSet = cat(4,img_4D, imgaussfilt3(single(img_4D(:,:,:,1)),1), X, Y, Z);    % imgaussfilt .*focus

                T = imsegkmeans3(featureSet, K, 'NormalizeInput',true);   %-- KMEANS --%

            otherwise
                error("unknown strategy")
        end

        % if(spatial_smoothing)   % to test
        %     neighb = 1;
        %     [mixstats.klas, K] = mrf_smoothing(coord, mixstats.klas, neighb);
        % end

        fprintf('Elapsed time %f sec \n', toc);


    end
    %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               RESULTS                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    
    clr = jet(K+round(1/3*K));
    clr = clr(1:K,:);  % remove red color

    vars.slic_show = slic_show;
    vars.slic_inds_1 = slic_inds(1);
    vars.subj_slic = subj_slic;
    vars.len_sp = len_sp;
    vars.K = K;
    vars.clr = clr;
    vars.mnx = [lvl-wdw/2, lvl+wdw/2];
    vars.max_slic_subplot = max_slic_subplot;
    vars.rmin = rmin;
    vars.rmax = rmax;
    vars.cmin = cmin;
    vars.cmax = cmax;
    vars.tumor_contour_list = tumor_contour_list;
      
        
    %%  Build label maps  &  Compute similarity scores

    switch(strategy)

        case('imsegkmeans')

            [dice_array, jacc_array] = compute_simi_scores(T, logical(segm_vol_full(rmin:rmax,cmin:cmax,slic_min+slic_inds-1)), K);
            [max_dice, max_cl] = max(dice_array);
            vars.max_cl = max_cl;
            
            if plot_rab
                figure(fig_slic6);
                show_result_on_img_kmeans(T, img_4D, focus(rmin:rmax,cmin:cmax,slic_show), vars)
                suptitle(sprintf('6 %s slices / %d.   Best cluster: Dice = %0.3f, IoU = %0.3f', show_slices, length(slic_inds), max_dice, jacc_array(max_cl)))
            end
            
            
            if ~isempty(fn_save_pdf)
                if ~load_saved_mdl
                    rdn = num2str(randi(1000));
                    save([fn_save_pdf,'_',rdn,'__ROI',num2str(roi_radius),'_cl',num2str(K),'_imsegkmeans.mat'],'T')
                end
                if plot_rab
                    savefig(fig_slic6,[fn_save_pdf,'_',rdn,'__ROI',num2str(roi_radius),'_cl',num2str(K),'_imsegkmeans']);
                    save_pdf(fig_slic6, [fn_save_pdf,'_',rdn,'__ROI',num2str(roi_radius),'_cl',num2str(K),'_imsegkmeans.pdf']);
                end
            end
            

        otherwise
        %%
        % save result variables
        
        if strcmp(strategy,'kmeans')
            if ~isempty(fn_save_pdf) && ~load_saved_mdl
                rdn = num2str(randi(1000));
                save([fn_save_pdf,'_',rdn,'__ROI',num2str(roi_radius),'_cl',num2str(K),'_Kmeans.mat'],'mixstats')
            end
        else
            if ~isempty(fn_save_pdf)  &&  ~load_saved_mdl
                rdn = num2str(randi(1000));
                save([fn_save_pdf,'_',rdn,'__ROI',num2str(roi_radius),'_cl',num2str(K),'_mixmodel.mat'],'mixmodel')
                save([fn_save_pdf,'_',rdn,'__ROI',num2str(roi_radius),'_cl',num2str(K),'_mixstats.mat'],'mixstats')
            end
        end


%         % Build tumor class label map
%         lin_tum = find(gr_truth);  % tumor;
%         klas_tum = ismember(sub2ind(size(segm_vol_full(:,:,slic_inds)), coord(:,1), coord(:,2), coord(:,3)-slic_inds(1)+1),lin_tum);
           
        % Reconstruct results label map as 3D volume
        reco_lbl = zeros(size(gr_truth));
        reco_lbl(lin_obj) = mixstats.klas;
        % Simi score
%         [dice_array, jacc_array] = compute_simi_scores(mixstats.klas, klas_tum, K);
%         [max_sim, max_cl] = max(dice_array);
        match_klas = match_tumor_clusters(klas_tum,mixstats.klas,K);
        [max_dice, max_jacc] = compute_best_simi_score(mixstats.klas, klas_tum, match_klas);
        
        % Apply median filter
        if plot_filter3
            fsz1 = 3; % filter size: odd nb
            reco_lbl_filt1 = medfilt3(reco_lbl, [fsz1 fsz1 fsz1]);
            % reco_lbl_filt = modefilt(reco_lbl, [fsz fsz fsz]);  % when using Matlab from 2020a
            klas_filt1 = reco_lbl_filt1(lin_obj);   % Reconstruct filtered label map as column vector
            % Simi score
%             [dice_array1, jacc_array1] = compute_simi_scores(klas_filt1, klas_tum, K);
%             [max_sim1, max_cl1] = max(dice_array1);
            match_klas_filt1 = match_tumor_clusters(klas_tum,klas_filt1,K);
            [max_dice_filt1, max_jacc_filt1] = compute_best_simi_score(klas_filt1, klas_tum, match_klas_filt1);
        end

        % Apply other median filter
        if plot_filter5
            fsz2 = 5; % filter size: odd nb
            reco_lbl_filt2 = medfilt3(reco_lbl, [fsz2 fsz2 fsz2]);
            % reco_lbl_filt = modefilt(reco_lbl, [fsz2 fsz2 fsz2]);  % when using Matlab from 2020a
            klas_filt2 = reco_lbl_filt2(lin_obj);   % Reconstruct filtered label map as column vector
            % Simi score
%             [dice_array2, jacc_array2] = compute_simi_scores(klas_filt2, klas_tum, K);
%             [max_sim2, max_cl2] = max(dice_array2);
            match_klas_filt2 = match_tumor_clusters(klas_tum,klas_filt2,K);
            [max_dice_filt2, max_jacc_filt2] = compute_best_simi_score(klas_filt2, klas_tum, match_klas_filt2);
        end
        
        
        %% Plot results on original image
        
        if plot_initResults
            figure(fig_slic)
%             vars.max_cl = max_cl;
            vars.match_klas = match_klas;
            show_result_on_img(reco_lbl, vars);
            suptitle(sprintf('6 %s slices / %d.   %d red clusters: Dice = %0.3f, IoU = %0.3f', show_slices, length(slic_inds), length(match_klas), max_dice, max_jacc))
        end
        
        % % With majority filter
        if plot_filter3
            figure(fig_slic2)
%             vars.max_cl = max_cl1;
            vars.match_klas = match_klas_filt1;
            show_result_on_img(reco_lbl_filt1, vars);
            suptitle(sprintf('6 %s slices / %d.   With %d*%d*%d filter  -  %d red clusters: Dice = %0.3f, IoU = %0.3f', show_slices, length(slic_inds), fsz1, fsz1, fsz1, length(match_klas), max_dice_filt1, max_jacc_filt1))
        end
        
        % % With majority filter 2
        if plot_filter5
            figure(fig_slic3)
%             vars.max_cl = max_cl2;
            vars.match_klas = match_klas_filt2;
            show_result_on_img(reco_lbl_filt2, vars);
            suptitle(sprintf('6 %s slices / %d.   With %d*%d*%d filter  -  %d red clusters: Dice = %0.3f, IoU = %0.3f', show_slices, length(slic_inds), fsz2, fsz2, fsz2, length(match_klas), max_dice_filt2, max_jacc_filt2))
        end
        
        if plot_rab
            figure(fig_slic6)
%             vars.max_cl = max_cl;
            vars.match_klas = match_klas_filt2;
            show_result_on_img(reco_lbl_filt2, vars);
            suptitle(sprintf('6 %s slices / %d.   With %d*%d*%d filter  -  %d red clusters: Dice = %0.3f, IoU = %0.3f', show_slices, length(slic_inds), fsz2, fsz2, fsz2, length(match_klas), max_sim2, jacc_array2(max_cl2)))
        end
        
        
        %% Merge close clusters
        
        klas = mixstats.klas;
        
        merged_cl = [0 0];
        for k = 1:K
            for k2 = (k+1):K
                % functional distance (between Muk)
                f = sum((mixstats.Muk(:,k2)-mixstats.Muk(:,k)).^2).^0.5 < 0.9;
                % spatial distance (between coord)
                coord_k = coord(klas == k,:); coord_k2 = coord(klas == k2,:);
                mk = median(coord_k,1); mk2 = median(coord_k2,1);
                dist_cl2_mk = mean(sum((mk-coord_k2).^2,2).^0.5);
                dist_cl_mk2 = mean(sum((mk2-coord_k).^2,2).^0.5);
                d = (dist_cl2_mk < 18) && (dist_cl_mk2 < 18);  %20;  % 32.5;    % d = norm(mk-mk2) < 20;
                
                % Merge
                if ( f && d )
                    merged_cl = [merged_cl; [k, k2]];
                    klas( klas == k ) = k2;
                end
            end
        end
        mixstats.merged_cl = merged_cl(2:end,:);
        
        
        %% Build label maps  &  Compute similarity scores  &  Plot results
        
        % reconstruct label map
        reco_lbl = zeros(size(gr_truth));
        reco_lbl(lin_obj) = klas;
        if plot_initMerge
            % compute similarity scores
%             [dice_array, jacc_array] = compute_simi_scores(klas, klas_tum, K);
%             [max_dice, max_cl] = max(dice_array);
%             vars.max_cl = max_cl;
            match_klas = match_tumor_clusters(klas_tum,klas,K);
            [max_dice, max_jacc] = compute_best_simi_score(klas, klas_tum, match_klas);
            vars.match_klas = match_klas;
            % plot results
            figure(fig_slic4)
            show_result_on_img(reco_lbl, vars);
            suptitle(sprintf(sprintf('6 %s slices / %d.   Merged  -  %d red clusters: Dice = %0.3f, IoU = %0.3f', show_slices, length(slic_inds), length(match_klas), max_dice, max_jacc)))
        end
        
        % apply filter
        if plot_mergeFilter3
            fsz1 = 3; % filter size: odd nb
            reco_lbl_filt1 = medfilt3(reco_lbl, [fsz1 fsz1 fsz1]);
            klas_filt1 = reco_lbl_filt1(lin_obj);
            % compute similarity scores
%             [dice_array1, jacc_array1] = compute_simi_scores(klas_filt1, klas_tum, K);
%             [max_sim1, max_cl1] = max(dice_array1);
%             vars.max_cl = max_cl1;
            match_klas_filt1 = match_tumor_clusters(klas_tum,klas_filt1,K);
            [max_dice_filt1, max_jacc_filt1] = compute_best_simi_score(klas_filt1, klas_tum, match_klas_filt1);
            vars.match_klas = match_klas_filt1;
            % plot results
            figure(fig_slic5)
            show_result_on_img(reco_lbl_filt1, vars);
            suptitle(sprintf('6 %s slices / %d.   Merged  -  With %d*%d*%d filter  -  %d red clusters: Dice = %0.3f, IoU = %0.3f', show_slices, length(slic_inds), fsz1, fsz1, fsz1, length(match_klas_filt1), max_dice_filt1, max_jacc_filt1))
        end
        
        
        %% Plot decay curve results
        
        if ~load_saved_mdl && plot_clusterCurves
            
            Curves.decay_curves = decay_curves;
            [fig_normCurves, fig_decayCurves, fig_loglik] = show_SRM_results_new(Curves, mixmodel, mixstats, model, max_cl, clr);
            
        end
        

    
        %% Save image in pdf
        if ~isempty(fn_save_pdf)   % && ~load_saved_mdl

            %     rdn = num2str(randi(1000));  % already initialized in SRM MODEL FITTING SECTION

            if plot_initResults, savefig(fig_slic,[fn_save_pdf,'_',rdn,'__ROI',num2str(roi_radius),'_cl',num2str(K)]); end
            if plot_filter3, savefig(fig_slic2,[fn_save_pdf,'_',rdn,'__ROI',num2str(roi_radius),'_cl',num2str(K),'_filter',num2str(fsz1)]); end
            if plot_filter5, savefig(fig_slic3,[fn_save_pdf,'_',rdn,'__ROI',num2str(roi_radius),'_cl',num2str(K),'_filter',num2str(fsz2)]); end
            if plot_initMerge, savefig(fig_slic4,[fn_save_pdf,'_',rdn,'__ROI',num2str(roi_radius),'_cl',num2str(K),'_merged']); end
            if plot_mergeFilter3, savefig(fig_slic5,[fn_save_pdf,'_',rdn,'__ROI',num2str(roi_radius),'_cl',num2str(K),'_merged_filter',num2str(fsz1)]); end

            if plot_initResults, save_pdf(fig_slic, [fn_save_pdf,'_',rdn,'__ROI',num2str(roi_radius),'_cl',num2str(K),'.pdf']); end
            if plot_filter3, save_pdf(fig_slic2,[fn_save_pdf,'_',rdn,'__ROI',num2str(roi_radius),'_cl',num2str(K),'_filter',num2str(fsz1),'.pdf']); end
            if plot_filter5, save_pdf(fig_slic3,[fn_save_pdf,'_',rdn,'__ROI',num2str(roi_radius),'_cl',num2str(K),'_filter',num2str(fsz2),'.pdf']); end
            if plot_initMerge, save_pdf(fig_slic4,[fn_save_pdf,'_',rdn,'__ROI',num2str(roi_radius),'_cl',num2str(K),'_merged']); end
            if plot_mergeFilter3, save_pdf(fig_slic5,[fn_save_pdf,'_',rdn,'__ROI',num2str(roi_radius),'_cl',num2str(K),'_merged_filter',num2str(fsz1)]); end


            % Save more figures

            if plot_clusterCurves
                save_pdf(fig_normCurves,[fn_save_pdf,'_',rdn,'__ROI',num2str(roi_radius),'_cl',num2str(K),'_mdlCurveDetails.pdf'])
                save_pdf(fig_decayCurves,[fn_save_pdf,'_',rdn,'__ROI',num2str(roi_radius),'_cl',num2str(K),'_decayCurveDetails.pdf'])
                save_pdf(fig_loglik,[fn_save_pdf,'_',rdn,'__ROI',num2str(roi_radius),'_cl',num2str(K),'_convergence.pdf'])
            end


            % save updated mixstats with merged_cluster array
%             save([fn_save_pdf,'_',rdn,'__ROI',num2str(roi_radius),'_cl',num2str(K),'_mixstats.mat'],'mixstats')
            
            % save reduced variable to be smaller in size and able to transfer easily from GPU to local machine
            mixstats_red.Muk = mixstats.Muk;
            mixstats_red.klas = mixstats.klas;
            save([fn_save_pdf,'_',rdn,'__ROI',num2str(roi_radius),'_cl',num2str(K),'_mixstats_red.mat'],'mixstats_red')
        end
    
    end
%     pat = [pat, patient_nm, ' '];
%     dic = [dic, max_dice_filt1];
%     jac = [jac, max_jacc_filt1];
    
end  % end patient

end  % end K

% pat, dic, jac



