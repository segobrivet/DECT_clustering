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
% machine_type = 'Peter';
% machine_type = 'GPU';
% machine_type = 'GPU_clean';

results_folder_name = 'results_sft';
% results_folder_name = 'results_test_kmeans';
% results_folder_name = 'results_kmeans';
% results_folder_name = 'results_square_150';
% results_folder_name = 'results_organ3';
% results_folder_name = 'results_clean';
% results_folder_name = 'results_gmm';

%%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
strategy = 'FunSRM'; init_kmeans_dep = 1;
% strategy = 'ExtraFeat-GMM';
% strategy = 'kmeans';
% strategy = 'gmm';
% strategy = 'imsegkmeans';

mixingOption = 'softmax';
% mixingOption = 'gaussian';

% model = "PRM"; % Polynomial regression mixture
% model = "SRM"; % Spline regression mixture
model = "bSRM";% B-Spline regression mixture

% p = 3; % polynomial regression degree
% spline_order = 4; % 2: linear, 3: qudratic, 4: cubic, etc, for (b-)spline spatial regression mixture
Bspline_order = 4; % 2: linear, 3: qudratic, 4: cubic, etc, for (b-)spline spatial regression mixture
nknots = 10; % fixed number of internal knots, for (b-)spline spatial regression mixture


%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROI size
% roi_radius = 90;   % for a circle  % default is 90
roi_radius = 100;    % for a square  % default is 150

% tissue enhancement window
lvl = 150;  % 40  % 70
wdw = 700;  % 350  % 500


org_ids = '3';  % 1-2

take_which_slices = 6;  % if any num: take num tumor slices in algo; if 'all': take all available slices (max 20 if latest generated) % default is 6

%%%%%%%%%%%%%%%%%%%%%%%%%%% CHOOSE PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_slic_subplot = 6;   % max nb of subplots on the same figure

show_slices = 'middle';  % if num slices > num subplots, which slices are we plotting?
% show_slices = 'bottom';
% show_slices = 'top';

plot_clusterCurves = 0;

plot_initResults = 1;
plot_filter3 = 1;
plot_filter5 = 0;
plot_initMerge = 0;
plot_mergeFilter3 = 0;
plot_rab = 0;
plot_subfigure = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%% TO LOAD A SAVED MODEL %%%%%%%%%%%%%%%%%%%%%%%%%
load_saved_mdl = 0;
pat = [];
dic = zeros(100,1);
jac = zeros(100,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%% NB OF RESTART %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nbr_EM_runs = 1;

%% end of specifications


% % % if load_saved_mdl % % %
%% To read saved result variables in a folder and plot the images
% fls = dir(fullfile(results_folder_name,'/HNSCC9_*__ROI150_cl40_lbda0.075_mixstats_red.mat'));
% 
% for nm = 1:length(fls)
% %     imsegkmeanspath = fullfile(results_folder_name,fls(nm).name);
%     mixstatspath = fullfile(results_folder_name,fls(nm).name);
% %     mixmodelpath = [mixstatspath(1:end-12),'mixmodel.mat'];  % (1:end-17)
%     str = split(fls(nm).name,'_');
%     patient_name = str{1}; rdn = str{2};
%     roi_radius = str{4}; roi_radius = str2num(roi_radius(4:end));
%     K = str{5}; K = str2num(K(3:end));
%     lbda = str{6}; lbda = str2num(lbda(5:end));
%     
% for lambda = lbda


% % for visualising init MyKmeans
% % fls = dir(fullfile(results_folder_name,'/solInitKmeans*54.mat'));
% fls = dir('solInitKmeans*815.mat');
% for nm = 1:length(fls)
% %     mixstatspath = fullfile(results_folder_name,fls(nm).name);
%     mixstatspath = fullfile(fls(nm).name);
%     str = split(fls(nm).name,'_');
%     patient_name = 'HNSCC3';
%     roi_radius = 90;
%     K = str{2}; K = str2num(K(3:end));
%     rdn = str{3}; rdn = rdn(1:end-4);
% for lambda = 0.07


% % for visualising GMM
% fls = dir(fullfile(results_folder_name,'/HNSCC10_*__gmm_ROI100_cl150_mixstats_red.mat'));
% for nm = 1:length(fls)
%     mixstatspath = fullfile(results_folder_name,fls(nm).name);
%     mixmodelpath = [mixstatspath(1:end-16),'mixmodel.mat'];
%     str = split(fls(nm).name,'_');
%     patient_name = str{1}; rdn = str{2};
%     roi_radius = str{5}; roi_radius = str2num(roi_radius(4:end));
%     K = str{6}; K = str2num(K(3:end));

%%

for lambda = [0.075]  % default is 0.075
    
pat_ind = 0;
for patient_name = ["HNSCC5","HNSCC8","HNSCC9","HNSCC10"] %,"HNSCC3","HNSCC5","HNSCC8","HNSCC9","HNSCC10"]
    
% for patient_name = ["HNSCC2","HNSCC3","HNSCC5","HNSCC8","HNSCC9","HNSCC10",...
%         "HNSCC11","HNSCC12","HNSCC13","HNSCC15","HNSCC15A","HNSCC17","HNSCC17A","HNSCC18","HNSCC20",...
%         "HNSCC21","HNSCC22","HNSCC22A","HNSCC25","HNSCC26","HNSCC27","HNSCC29","HNSCC30",...
%         "HNSCC31A","HNSCC32","HNSCC33","HNSCC34","HNSCC35","HNSCC36","HNSCC37A","HNSCC38","HNSCC39",...
%         "HNSCC41","HNSCC42","HNSCC44","HNSCC44AM","HNSCC45","HNSCC46","HNSCC47","HNSCC48","HNSCC49",...
%         "HNSCC51","HNSCC52","HNSCC52AM","HNSCC53","HNSCC55","HNSCC56","HNSCC57",...
%         "HNSCC61A","HNSCC62","HNSCC63","HNSCC63A","HNSCC64A","HNSCC65A","HNSCC66","HNSCC67","HNSCC68","HNSCC69",...
%         "HNSCC70A","HNSCC71","HNSCC72A","HNSCC73","HNSCC74","HNSCC75","HNSCC76","HNSCC77","HNSCC78","HNSCC79","HNSCC80",...
%         "HNSCC81","HNSCC82","HNSCC83","HNSCC84","HNSCC85","HNSCC87","HNSCC88","HNSCC89","HNSCC90",...
%         "HNSCC91","HNSCC92","HNSCC95","HNSCC96","HNSCC97","HNSCC98",...
%         "HNSCC100","HNSCC101","HNSCC103","HNSCC105","HNSCC106","HNSCC108","HNSCC109"]
% 	% data with pb, taken out from all available patients above: "HNSCC1","HNSCC10A","HNSCC60","HNSCC102"
        
    pat_ind = pat_ind + 1;
        
    
nb_K_ind = 0;
for K = [150]      % default for FunSRM is K=40 with square ROI of size 150   % default for imsegkmeans: 150
    
    nb_K_ind = nb_K_ind + 1;
    
    
    %% Data (spectral image)

    close all; 
    
    patient_nm = char(patient_name); disp(patient_nm);
%     fn_save_pdf = '';  % if don't want to save fig as pdf
    fn_save_pdf = fullfile(results_folder_name, patient_nm);  % to save fig as pdf
    
    %%% Build dataset for a specific patient
    
    res = build_patient_data(machine_type, patient_nm, take_which_slices, org_ids, roi_radius, max_slic_subplot, show_slices);
    
    Y = res.Y; decay_curves = res.decay_curves; gr_truth = res.gr_truth;
    klas_tum = res.klas_tum; lin_obj = res.lin_obj; subj_slic = res.subj_slic; focus = res.focus;
    slic_show = res.slic_show; slic_min_idx = res.slic_min_idx; slic_min = res.slic_min; slic_inds = res.slic_inds;
    rmin = res.rmin; rmax = res.rmax; cmin = res.cmin; cmax = res.cmax; coord = res.coord;
    
    [n, m] = size(Y);
    V = coord; V(:,3) = coord(:,3)*2;   % spacing in Z is 1.25mm, spacing in X and in Y is 0.61mm, so Z is 2*bigger than X-Y.
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
                
            case('gmm')
                load(mixstatspath);
                load(mixmodelpath);
                mixstats.klas = mixstats_red.klas;
                mixstats.Muk = mixmodel.mu;
                
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

                
            case('gmm')
                mixmodel = fitgmdist(Curves.ordinates,K,'RegularizationValue',0.1,...
                                     'Options',statset('MaxIter',1500,'TolFun',1e-6));
                [mixstats.klas, mixstats.negloglik, mixstats.posterior, mixstats.logpdf, mixstats.mahadist] = ...
                                                                                  cluster(mixmodel,Curves.ordinates);
                
                
            otherwise
                error("unknown strategy")
        end

        % if(spatial_smoothing)   % to test
        %     neighb = 1;
        %     [mixstats.klas, K] = mrf_smoothing(coord, mixstats.klas, neighb);
        % end

        elaps_time_clust = toc;
        fprintf('Elapsed time for clustering algo: %f sec \n', elaps_time_clust);


    end
    %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               RESULTS                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    
    clr = jet(K+round(1/3*K));
    clr = clr(1:K,:);  % this removes red color
    
    nb_subplot = min(max_slic_subplot,length(slic_show));
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
        elseif strcmp(strategy,'gmm')
            if ~isempty(fn_save_pdf) && ~load_saved_mdl
                rdn = num2str(randi(1000));
                save([fn_save_pdf,'_',rdn,'__ROI',num2str(roi_radius),'_cl',num2str(K),'_GMM.mat'],'mixstats')
            end
        else
            if ~isempty(fn_save_pdf)  &&  ~load_saved_mdl
                rdn = num2str(randi(1000));
                save([fn_save_pdf,'_',rdn,'__ROI',num2str(roi_radius),'_cl',num2str(K),'_lbda',num2str(lambda),'_mixmodel.mat'],'mixmodel')
                save([fn_save_pdf,'_',rdn,'__ROI',num2str(roi_radius),'_cl',num2str(K),'_lbda',num2str(lambda),'_mixstats.mat'],'mixstats')
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
%         match_klas = match_tumor_clusters(klas_tum,mixstats.klas,K);
%         [max_dice, max_jacc] = compute_best_simi_score(mixstats.klas, klas_tum, match_klas);
        [max_dice, max_jacc, match_klas] = compute_best_dicjac(klas_tum, mixstats.klas, K);
        
        % Apply median filter
        if plot_filter3
            fsz1 = 3; % filter size: odd nb
            reco_lbl_filt1 = medfilt3(reco_lbl, [fsz1 fsz1 fsz1]);
            % reco_lbl_filt = modefilt(reco_lbl, [fsz fsz fsz]);  % when using Matlab from 2020a
            klas_filt1 = reco_lbl_filt1(lin_obj);   % Reconstruct filtered label map as column vector
            % Simi score
%             [dice_array1, jacc_array1] = compute_simi_scores(klas_filt1, klas_tum, K);
%             [max_sim1, max_cl1] = max(dice_array1);
%             match_klas_filt1 = match_tumor_clusters(klas_tum,klas_filt1,K);
%             [max_dice_filt1, max_jacc_filt1] = compute_best_simi_score(klas_filt1, klas_tum, match_klas_filt1);
            [max_dice_filt1, max_jacc_filt1, match_klas_filt1] = compute_best_dicjac(klas_tum, klas_filt1, K);
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
%             match_klas_filt2 = match_tumor_clusters(klas_tum,klas_filt2,K);
%             [max_dice_filt2, max_jacc_filt2] = compute_best_simi_score(klas_filt2, klas_tum, match_klas_filt2);
            [max_dice_filt2, max_jacc_filt2, match_klas_filt2] = compute_best_dicjac(klas_tum, klas_filt2, K);
        end
        
        
        %% Plot results on original image
                
        if plot_initResults
            figure(fig_slic)
%             vars.max_cl = max_cl;
            vars.match_klas = match_klas;
            show_result_on_img(reco_lbl, vars);
            sgtitle(sprintf('%d slices / %d.   %d red clusters: Dice = %0.3f, IoU = %0.3f', nb_subplot, length(slic_inds), length(match_klas), max_dice, max_jacc),'FontSize',14,'FontWeight','bold');
        end
        
        % % With majority filter
        if plot_filter3
            figure(fig_slic2)
%             vars.max_cl = max_cl1;
            vars.match_klas = match_klas_filt1;
            show_result_on_img(reco_lbl_filt1, vars);
            sgtitle(sprintf('%d slices / %d.   With %d*%d*%d filter  -  %d red clusters: Dice = %0.3f, IoU = %0.3f', nb_subplot, length(slic_inds), fsz1, fsz1, fsz1, length(match_klas), max_dice_filt1, max_jacc_filt1),'FontSize',14,'FontWeight','bold');
        end
        
        % % With majority filter 2
        if plot_filter5
            figure(fig_slic3)
%             vars.max_cl = max_cl2;
            vars.match_klas = match_klas_filt2;
            show_result_on_img(reco_lbl_filt2, vars);
            sgtitle(sprintf('%d slices / %d.   With %d*%d*%d filter  -  %d red clusters: Dice = %0.3f, IoU = %0.3f', nb_subplot, length(slic_inds), fsz2, fsz2, fsz2, length(match_klas), max_dice_filt2, max_jacc_filt2),'FontSize',14,'FontWeight','bold');
        end
        
        if plot_rab
            figure(fig_slic6)
%             vars.max_cl = max_cl;
            vars.match_klas = match_klas_filt2;
            show_result_on_img(reco_lbl_filt2, vars);
            sgtitle(sprintf('%d slices / %d.   With %d*%d*%d filter  -  %d red clusters: Dice = %0.3f, IoU = %0.3f', nb_subplot, length(slic_inds), fsz2, fsz2, fsz2, length(match_klas), max_sim2, jacc_array2(max_cl2)),'FontSize',14,'FontWeight','bold');
        end
        
        
        %% Merge close clusters
        
%         klas = mixstats.klas;
%         
%         merged_cl = [0 0];
%         for k = 1:K
%             for k2 = (k+1):K
%                 % functional distance (between Muk)
%                 f = sum((mixstats.Muk(:,k2)-mixstats.Muk(:,k)).^2).^0.5 < 0.9;
%                 % spatial distance (between coord)
%                 coord_k = coord(klas == k,:); coord_k2 = coord(klas == k2,:);
%                 mk = median(coord_k,1); mk2 = median(coord_k2,1);
%                 dist_cl2_mk = mean(sum((mk-coord_k2).^2,2).^0.5);
%                 dist_cl_mk2 = mean(sum((mk2-coord_k).^2,2).^0.5);
%                 d = (dist_cl2_mk < 18) && (dist_cl_mk2 < 18);  %20;  % d = norm(mk-mk2) < 20;  %%% FIXED PARAM %%%
%                 
%                 % Merge
%                 if ( f && d )
%                     merged_cl = [merged_cl; [k, k2]];
%                     klas( klas == k ) = k2;
%                 end
%             end
%         end
%         mixstats.merged_cl = merged_cl(2:end,:);
        
        
        %% Build label maps  &  Compute similarity scores  &  Plot results
        
        %%% NOT UP TO DATE with 'compute_best_dicjac' function
        
%         % reconstruct label map
%         reco_lbl = zeros(size(gr_truth));
%         reco_lbl(lin_obj) = klas;
%         if plot_initMerge
%             % compute similarity scores
% %             [dice_array, jacc_array] = compute_simi_scores(klas, klas_tum, K);
% %             [max_dice, max_cl] = max(dice_array);
% %             vars.max_cl = max_cl;
%             match_klas = match_tumor_clusters(klas_tum,klas,K);
%             [max_dice, max_jacc] = compute_best_simi_score(klas, klas_tum, match_klas);
%             vars.match_klas = match_klas;
%             % plot results
%             figure(fig_slic4)
%             show_result_on_img(reco_lbl, vars);
%             suptitle(sprintf(sprintf('%d slices / %d.   Merged  -  %d red clusters: Dice = %0.3f, IoU = %0.3f', nb_subplot, length(slic_inds), length(match_klas), max_dice, max_jacc)))
%         end
%         
%         % apply filter
%         if plot_mergeFilter3
%             fsz1 = 3; % filter size: odd nb
%             reco_lbl_filt1 = medfilt3(reco_lbl, [fsz1 fsz1 fsz1]);
%             klas_filt1 = reco_lbl_filt1(lin_obj);
%             % compute similarity scores
% %             [dice_array1, jacc_array1] = compute_simi_scores(klas_filt1, klas_tum, K);
% %             [max_sim1, max_cl1] = max(dice_array1);
% %             vars.max_cl = max_cl1;
%             match_klas_filt1 = match_tumor_clusters(klas_tum,klas_filt1,K);
%             [max_dice_filt1, max_jacc_filt1] = compute_best_simi_score(klas_filt1, klas_tum, match_klas_filt1);
%             vars.match_klas = match_klas_filt1;
%             % plot results
%             figure(fig_slic5)
%             show_result_on_img(reco_lbl_filt1, vars);
%             sgtitle(sprintf('%d slices / %d.   Merged  -  With %d*%d*%d filter  -  %d red clusters: Dice = %0.3f, IoU = %0.3f', nb_subplot, length(slic_inds), fsz1, fsz1, fsz1, length(match_klas_filt1), max_dice_filt1, max_jacc_filt1),'FontSize',14,'FontWeight','bold')
%         end
        
        
        %% Plot decay curve results
        
        if ~load_saved_mdl && plot_clusterCurves
            
            Curves.decay_curves = decay_curves;
            [fig_normCurves, fig_decayCurves, fig_loglik] = show_SRM_results_new(Curves, mixmodel, mixstats, model, match_klas, clr);
            
            save_pdf(gcf,[fn_save_pdf,'__decay_curves.pdf']);
            
%             save_pdf(gcf,[fn_save_pdf,'_',rdn,'__ROI',num2str(roi_radius),'_cl',num2str(K),'_curves_4.pdf']);
        end
        
        if plot_subfigure
            i = 3;
            
            fig_slicOrig = figure('units','pixels','outerposition',[0 0 1500 1500]); 
            imshow(slics(:,:,slic_show(i)),[])
            hold on
%             plot_tumor_contour(tumor_contour_list{i}, [rmin, cmin], [0,0.5,1]);
            save_pdf(fig_slicOrig,[fn_save_pdf,'__ROI',num2str(roi_radius),'_onlySlic',num2str(i),'.pdf']);

        
            fig_subfig = figure('units','pixels','outerposition',[0 0 900 900]); 
            vars.match_klas = match_klas_filt1;
            show_result_one_subfigure(reco_lbl_filt1, vars, i);
%             save_pdf(fig_subfig,[fn_save_pdf,'_',rdn,'__ROI',num2str(roi_radius),'_cl',num2str(K),...
%                                 '_lbda',num2str(lambda),'_filter',num2str(fsz1),'onlySlic',num2str(i),'.pdf']);
        end

        
        %% Save image in pdf
        if ~isempty(fn_save_pdf)   % && ~load_saved_mdl

            if strcmp(strategy,'gmm')
                mixstats_red.klas = mixstats.klas;
                if exist('elaps_time_clust','var'), mixstats_red.elaps_time = elaps_time_clust; else, mixstats_red.elaps_time = NaN; end
                mixstats_red.dice = max_dice;
                if exist('max_dice_filt1','var'), mixstats_red.dice_mrg_flt = max_dice_filt1; else, mixstats_red.dice_mrg_flt = NaN; end
                save([fn_save_pdf,'_',rdn,'__gmm_ROI',num2str(roi_radius),'_cl',num2str(K),'_mixstats_red.mat'],'mixstats_red')
                save([fn_save_pdf,'_',rdn,'__gmm_ROI',num2str(roi_radius),'_cl',num2str(K),'_mixstats.mat'],'mixstats')
                save([fn_save_pdf,'_',rdn,'__gmm_ROI',num2str(roi_radius),'_cl',num2str(K),'_mixmodel.mat'],'mixmodel')
            
                if plot_initResults, save_pdf(fig_slic, [fn_save_pdf,'_',rdn,'__gmm_ROI',num2str(roi_radius),'_cl',num2str(K),'.pdf']); end
            
            else
                
            
            %     rdn = num2str(randi(1000));  % already initialized in SRM MODEL FITTING SECTION

            % Save mixtats with more information
            
            % save updated mixstats with merged_cluster array
%             save([fn_save_pdf,'_',rdn,'__ROI',num2str(roi_radius),'_cl',num2str(K),'_mixstats.mat'],'mixstats')
            
            % save reduced variable to be smaller in size and able to transfer easily from GPU to local machine
            mixstats_red.Muk = mixstats.Muk;
            mixstats_red.klas = mixstats.klas;
            mixstats_red.lambda = lambda;
            if exist('elaps_time_clust','var'), mixstats_red.elaps_time = elaps_time_clust; else, mixstats_red.elaps_time = NaN; end
            mixstats_red.loglik = mixstats.loglik;
            mixstats_red.dice = max_dice;
            if exist('max_dice_filt1','var'), mixstats_red.dice_mrg_flt = max_dice_filt1; else, mixstats_red.dice_mrg_flt = NaN; end
            save([fn_save_pdf,'_',rdn,'__ROI',num2str(roi_radius),'_cl',num2str(K),'_lbda',num2str(lambda),'_mixstats_red.mat'],'mixstats_red')
            
            
            % Save images
            
            if plot_initResults, savefig(fig_slic,[fn_save_pdf,'_',rdn,'__ROI',num2str(roi_radius),'_cl',num2str(K),'_lbda',num2str(lambda),'.fig']); end
            if plot_filter3, savefig(fig_slic2,[fn_save_pdf,'_',rdn,'__ROI',num2str(roi_radius),'_cl',num2str(K),'_lbda',num2str(lambda),'_filter',num2str(fsz1),'.fig']); end
            if plot_filter5, savefig(fig_slic3,[fn_save_pdf,'_',rdn,'__ROI',num2str(roi_radius),'_cl',num2str(K),'_lbda',num2str(lambda),'_filter',num2str(fsz2),'.fig']); end
            if plot_initMerge, savefig(fig_slic4,[fn_save_pdf,'_',rdn,'__ROI',num2str(roi_radius),'_cl',num2str(K),'_lbda',num2str(lambda),'_merged','.fig']); end
            if plot_mergeFilter3, savefig(fig_slic5,[fn_save_pdf,'_',rdn,'__ROI',num2str(roi_radius),'_cl',num2str(K),'_lbda',num2str(lambda),'_merged_filter',num2str(fsz1),'.fig']); end

            if plot_initResults, save_pdf(fig_slic, [fn_save_pdf,'_',rdn,'__ROI',num2str(roi_radius),'_cl',num2str(K),'_lbda',num2str(lambda),'.pdf']); end
            if plot_filter3, save_pdf(fig_slic2,[fn_save_pdf,'_',rdn,'__ROI',num2str(roi_radius),'_cl',num2str(K),'_lbda',num2str(lambda),'_filter',num2str(fsz1),'.pdf']); end
            if plot_filter5, save_pdf(fig_slic3,[fn_save_pdf,'_',rdn,'__ROI',num2str(roi_radius),'_cl',num2str(K),'_lbda',num2str(lambda),'_filter',num2str(fsz2),'.pdf']); end
            if plot_initMerge, save_pdf(fig_slic4,[fn_save_pdf,'_',rdn,'__ROI',num2str(roi_radius),'_cl',num2str(K),'_lbda',num2str(lambda),'_merged','.pdf']); end
            if plot_mergeFilter3, save_pdf(fig_slic5,[fn_save_pdf,'_',rdn,'__ROI',num2str(roi_radius),'_cl',num2str(K),'_lbda',num2str(lambda),'_merged_filter',num2str(fsz1),'.pdf']); end


            % Save more figures

            if plot_clusterCurves
                save_pdf(fig_normCurves,[fn_save_pdf,'_',rdn,'__ROI',num2str(roi_radius),'_cl',num2str(K),'_mdlCurveDetails.pdf'])
                save_pdf(fig_decayCurves,[fn_save_pdf,'_',rdn,'__ROI',num2str(roi_radius),'_cl',num2str(K),'_decayCurveDetails.pdf'])
                save_pdf(fig_loglik,[fn_save_pdf,'_',rdn,'__ROI',num2str(roi_radius),'_cl',num2str(K),'_convergence.pdf'])
            end

            end
        end
    
    end
    
    dic(pat_ind,nb_K_ind) = max_dice;
    jac(pat_ind,nb_K_ind) = max_jacc;

end  % end K
    
    pat = [pat, patient_nm, '; '];

end  % end patient

end  % end lambda

save([results_folder_name,'/all_dic.mat'],'dic');
save([results_folder_name,'/all_jac.mat'],'jac');
save([results_folder_name,'/all_pat.mat'],'pat');
% pat, dic, jac



