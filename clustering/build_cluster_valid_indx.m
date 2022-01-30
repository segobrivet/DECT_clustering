

%%% validity indices for clustering

addpath('utils');

clear, clc

%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

machine_type = 'sego_clean';

results_folder_name = 'results-SgMFR-Bspl';
% results_folder_name = 'results-SgMFR-poly';
% results_folder_name = 'results-SgMVFR-Bspl';
% results_folder_name = 'results-SsMFR';
% results_folder_name = 'results-GMM';
% results_folder_name = 'results-kmeans';
% results_folder_name = 'results-SelSearch';

strategy = 'ours';
% strategy = 'gmm';
% strategy = 'kmeans';
% strategy = 'selsearch';

% tissue enhancement window
lvl = 150;  % 40  % 70
wdw = 700;  % 350  % 500

plot_Stats = 0;
plot_Results = 0;
nb_slices = 6;  % % if any num: take num slices; if 'all': take all available slices (max 20 if latest generated)
max_slic_subplot = 6;
show_slices = 'middle';

pat = [];
nb_K = zeros(1,91);
dic = zeros(8,91);
jac = zeros(8,91);

db_spat = zeros(8,91);
ch_spat = zeros(8,91);
kl_spat = zeros(8,91);

db_spec = zeros(8,91);
ch_spec = zeros(8,91);
kl_spec = zeros(8,91);

db_spatec = zeros(9,91);

final_resuls_table = zeros(6,91);  % 1.dice 2.DBspat 3.DBspec 4.DBt_spat 5. DBt_spec, 6.time

%%%%%%%%%%%%%%%%%%%%%%%%%%% READ SAVED DATA %%%%%%%%%%%%%%%%%%%%%%%%%
roi_radius = 150;

crit_idx = 0;

patient_names = ["HNSCC76"]; %,"HNSCC5","HNSCC8","HNSCC9","HNSCC10"]
    
% patient_names = ["HNSCC2","HNSCC3","HNSCC5","HNSCC8","HNSCC9","HNSCC10",...
%         "HNSCC11","HNSCC12","HNSCC13","HNSCC15","HNSCC15A","HNSCC17","HNSCC17A","HNSCC18","HNSCC20",...
%         "HNSCC21","HNSCC22","HNSCC22A","HNSCC25","HNSCC26","HNSCC27","HNSCC29","HNSCC30",...
%         "HNSCC31A","HNSCC32","HNSCC33","HNSCC34","HNSCC35","HNSCC36","HNSCC37A","HNSCC38","HNSCC39",...
%         "HNSCC41","HNSCC42","HNSCC44","HNSCC44AM","HNSCC45","HNSCC46","HNSCC47","HNSCC48","HNSCC49",...
%         "HNSCC51","HNSCC52","HNSCC52AM","HNSCC53","HNSCC55","HNSCC56","HNSCC57",...
%         "HNSCC61A","HNSCC62","HNSCC63","HNSCC63A","HNSCC64A","HNSCC65A","HNSCC66","HNSCC67","HNSCC68","HNSCC69",...
%         "HNSCC70A","HNSCC71","HNSCC72A","HNSCC73","HNSCC74","HNSCC75","HNSCC76","HNSCC77","HNSCC78","HNSCC79","HNSCC80",...
%         "HNSCC81","HNSCC82","HNSCC83","HNSCC84","HNSCC85","HNSCC87","HNSCC88","HNSCC89","HNSCC90",...
%         "HNSCC91","HNSCC92","HNSCC95","HNSCC96","HNSCC97","HNSCC98",...
%         "HNSCC100","HNSCC101","HNSCC103","HNSCC105","HNSCC106","HNSCC108","HNSCC109"];
% 		% "HNSCC1","HNSCC10A","HNSCC60","HNSCC102"]
for patient_name = patient_names

    for K = 40 % [30, 35, 40, 50, 60, 70]

        crit_idx = crit_idx + 1;
    
    %% Load subject scan
    
    patient_nm = char(patient_name); disp(patient_nm);
    % fn_save_pdf = '';  % if don't want to save fig as pdf
    fn_save_pdf = fullfile(results_folder_name, patient_nm);  % to save fig as pdf
    
    res = build_patient_data(machine_type, patient_nm, nb_slices, roi_radius, max_slic_subplot, show_slices);
    
    Y = res.Y; decay_curves = res.decay_curves; gr_truth = res.gr_truth;
    klas_tum = res.klas_tum; lin_obj = res.lin_obj; subj_slic = res.subj_slic; focus = res.focus;
    slic_show = res.slic_show; slic_min_idx = res.slic_min_idx; slic_min = res.slic_min; slic_inds = res.slic_inds;
    rmin = res.rmin; rmax = res.rmax; cmin = res.cmin; cmax = res.cmax; coord = res.coord;
    
    [n, m] = size(Y);
    V = coord; V(:,3) = coord(:,3)*2;   % spacing in Z is 1.25mm, spacing in X and in Y is 0.61mm, so Z is 2*bigger than X-Y.
    T = linspace(0, 1, 21);  % T = linspace(40,140,21); %1:21;
    
    
    %% Pick param 
    
%     ticklbls = {'l=0.055','l=0.070','l=0.075','l=0.080','l=0.100'};
%     ticklbls = {'l=0.075','l=0.090','l=0.095','l=0.100','l=0.150','l=0.200'};
    ticklbls = {'SsFMR'};
%     ticklbls = {'K = 40','K = 150'};

    lbda_idx = 1;  %!%
    lambda = 0.075;
%     for lambda =  [0.055, 0.07, 0.075 0.08, 0.1] %[0.075, 0.09, 0.095, 0.1, 0.15, 0.2]  % [0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.09, 0.2]

%         lbda_idx = lbda_idx + 1;
        
        %% To read saved result variables in a folder
%         fls = dir(fullfile(results_folder_name,[patient_nm,'_*__ROI',num2str(roi_radius),'_cl',num2str(K),'_lbda',num2str(lambda),'_mixstats_red.mat']));
%         fls = dir(fullfile(results_folder_name,[patient_nm,'_*_ROI',num2str(roi_radius),'_cl',num2str(K),'*_mixstats_red.mat']));
%         fls = dir(fullfile(results_folder_name,'HNSCC15A_815__ROI150_cl40_lbda0.08_mixstats_red.mat'));
        fls = dir(fullfile(results_folder_name,[patient_nm,'_*_mixstats_red.mat']));
        if isempty(fls)
            db_spatec(:,crit_idx) = NaN;
            final_resuls_table(:,crit_idx) = NaN;
            continue;
        end
        
        nm = 1; %for nm = 1:length(fls)
            
            mixstatspath = fullfile(results_folder_name,fls(nm).name);
            
            %% Load saved model
            
            load(mixstatspath);

            if strcmp(strategy,'imsegkmeans')

                % % Segmentation score
                [max_dice, max_jacc, match_klas] = compute_best_dicjac(logical(gr_truth(rmin:rmax,cmin:cmax,:)), mixstats_red.T, K);
                
                % % Merge tumorous clusters
                lbl = mixstats_red.T;
                for kt=2:length(match_klas)
                    it = find(mixstats_red.T == match_klas(kt));
                    lbl(it) = match_klas(1);
                end
                % % remove clusters outside of focus (air data)
                lbl(~focus(rmin:rmax,cmin:cmax,:)) = 0;
                % % transpose to vector label corresponding to coordinates in V
                labels_mat = zeros(size(res.gr_truth));
                labels_mat(rmin:rmax,cmin:cmax,:) = lbl;
                elem_coord = res.coord;
                elem_coord(:,3) = elem_coord(:,3)-slic_inds(1)+1;
                labels = zeros(size(elem_coord,1),1);
                for i=1:size(elem_coord,1)
                	labels(i) = labels_mat(elem_coord(i,1), elem_coord(i,2), elem_coord(i,3));
                end

            elseif strcmp(strategy,'selsearch')
                
                str = split(fls(nm).name,'_');
                z_slic = str{end-2}; z_slic = str2double(z_slic(end));
                
                K = length(unique(mixstats_red.T));
                nb_K(crit_idx) = K;
                % % Segmentation score
                [max_dice, max_jacc, match_klas] = compute_best_dicjac(logical(gr_truth(rmin:rmax,cmin:cmax,z_slic)), mixstats_red.T, K);
                
                % % Merge tumorous clusters
                lbl = mixstats_red.T;
                for kt=2:length(match_klas)
                    it = find(mixstats_red.T == match_klas(kt));
                    lbl(it) = match_klas(1);
                end
                % % remove clusters outside of focus (air data)
                lbl(~focus(rmin:rmax,cmin:cmax,3)) = 0;
                % % transpose to vector label corresponding to coordinates in V
                labels_mat = zeros(size(res.gr_truth));
                labels_mat(rmin:rmax,cmin:cmax,3) = lbl;
                elem_coord = res.coord;
                elem_coord(:,3) = elem_coord(:,3)-slic_inds(1)+1;
                labels = zeros(size(elem_coord,1),1);
                for i=1:size(elem_coord,1)
                	labels(i) = labels_mat(elem_coord(i,1), elem_coord(i,2), elem_coord(i,3));
                end


            else
                if ~strcmp(strategy,'gmm')
                    mean_curves = mixstats_red.Muk;
                end
                
                % % Reconstruct results label map as 3D volume
                reco_lbl = zeros(size(gr_truth));
                reco_lbl(lin_obj) = mixstats_red.klas;
                % Apply median filter
                fsz1 = 3; % filter size: odd nb
                reco_lbl_filt1 = medfilt3(reco_lbl, [fsz1 fsz1 fsz1]);
                % reco_lbl_filt = modefilt(reco_lbl, [fsz fsz fsz]);  % when using Matlab from 2020a
                klas_filt1 = reco_lbl_filt1(lin_obj);   % Reconstruct filtered label map as column vector
            
                
                K = length(unique(mixstats_red.klas));
                nb_K(crit_idx) = K;
                
                % % Segmentation score
                [max_dice, max_jacc, match_klas] = compute_best_dicjac(klas_tum, klas_filt1, K);
                
                % % Merge tumorous clusters
                labels = klas_filt1;
                for kt=2:length(match_klas)
                    it = find(klas_filt1 == match_klas(kt));
                    labels(it) = match_klas(1);
                end
            
            end
            
            
            dic(lbda_idx,crit_idx) = max_dice;
            jac(lbda_idx,crit_idx) = max_jacc;
            
            final_resuls_table(1,crit_idx) = max_dice;
            
            %% Compute validity indices
            
            % % Spatial
            data = V;
            
            [indx,ssw,sw,sb] = valid_clusterIndex(data,labels,[]);
            
            db_spat(lbda_idx,crit_idx) = indx(2);
            ch_spat(lbda_idx,crit_idx) = indx(3);
            kl_spat(lbda_idx,crit_idx) = indx(4);
            
            db_spat_list = valid_DB_index(V,labels, match_klas(1));
            
            db_spatec(1,crit_idx) = nanmedian(db_spat_list);
            db_spatec(2,crit_idx) = nanmean(db_spat_list);
            db_spatec(3,crit_idx) = nanmax(db_spat_list);
            db_spatec(4,crit_idx) = nanstd(db_spat_list);
            
            
            final_resuls_table(2,crit_idx) = indx(2);
            final_resuls_table(4,crit_idx) = nanmax(db_spat_list);
            
            % % Curves
            data = Y;
            
            if ~strcmp(strategy,'ours') || size(mean_curves,1) ~= 21
                [indx,ssw,sw,sb] = valid_clusterIndex(data,labels,[]);
            else
                [indx,ssw,sw,sb] = valid_clusterIndex(data,labels,mean_curves');
            end

            db_spec(lbda_idx,crit_idx) = indx(2);
            ch_spec(lbda_idx,crit_idx) = indx(3);
            kl_spec(lbda_idx,crit_idx) = indx(4);
            
            
            db_spec_list = valid_DB_index(Y,labels, match_klas(1));
            
            db_spatec(5,crit_idx) = nanmedian(db_spec_list);
            db_spatec(6,crit_idx) = nanmean(db_spec_list);
            db_spatec(7,crit_idx) = nanmax(db_spec_list);
            db_spatec(8,crit_idx) = nanstd(db_spec_list);
            
            
            final_resuls_table(3,crit_idx) = indx(2);
            final_resuls_table(5,crit_idx) = nanmax(db_spec_list);
            
            
            % % Time
            final_resuls_table(6,crit_idx) = mixstats_red.elaps_time;
            
            
            mixstats_red.dice = final_resuls_table(1,crit_idx);
            mixstats_red.DB_spat = final_resuls_table(2,crit_idx);
            mixstats_red.DB_spec = final_resuls_table(3,crit_idx);
            mixstats_red.DBc_spat = final_resuls_table(4,crit_idx);
            mixstats_red.DBc_spec = final_resuls_table(5,crit_idx);
            
            save(mixstatspath, 'mixstats_red')
            
            %% Plot original image with tissue enhancement
            
            if plot_Results
                
                for i=1:length(slic_show)
                    tumor_contour_list{i} = bwboundaries(gr_truth(:,:,slic_show(i)));
                end
                
                clr = jet(K+round(1/3*K));
                clr = clr(1:K,:);  % this removes red color
                
                subj_slic(subj_slic<(lvl-(wdw/2))) = lvl-(wdw/2);
                subj_slic(subj_slic>(lvl+(wdw/2))) = lvl+(wdw/2);
                len_sp = min(ceil(length(slic_show)/2),ceil(max_slic_subplot/2));
                slics = subj_slic(rmin:rmax,cmin:cmax,:);
                
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
                
                fig_slic = figure('units','normalized','outerposition',[0 0 1 1]);
                for i=1:length(slic_show)
                    if i > max_slic_subplot, break; end
                    subplot(len_sp,4,2*i-1);
                    imshow(slics(:,:,slic_show(i)),[]), hold on, title("slice " + num2str(slic_inds(1)+slic_show(i)-1+slic_min-2)) % index showed on 3DSlicer
                    plot_tumor_contour(tumor_contour_list{i}, [rmin, cmin], [0,0.5,1]);
                end
                vars.match_klas = match_klas;
                show_result_on_img(reco_lbl_filt1, vars);
                sgtitle(sprintf('%d red clusters: Dice = %0.3f, DB = %0.3f / %0.3f, DB_c = %0.3f / %0.3f', length(vars.match_klas),...
                    max_dice, final_resuls_table(1,crit_idx), final_resuls_table(2,crit_idx), final_resuls_table(3,crit_idx), final_resuls_table(4,crit_idx)),'FontSize',14,'FontWeight','bold')
%                 save_pdf(gcf, [mixstatspath(1:end-17),'.pdf']);
            
            end
            
%         end
        
    end
    
end


% % % % % % % % % % % %  PLOT STATS  % % % % % % % % % % % % % % % % % % %
if plot_Stats
    n = length(ticklbls);
    figure;
    %%%
    subplot(3,1,1), hold on
    plot(db_spat(1:n,1),'b*--');
    plot(db_spec(1:n,1),'rx--');
    xticks(1:n); xticklabels(ticklbls);
    xlim([0.9, n+0.5])
    title('DB validity index')
    legend({'Spatial','Spectral'})
    %%%
    subplot(3,1,2), hold on
    plot(ch_spat(1:n,1),'b*--');
    plot(ch_spec(1:n,1),'rx--');
    xticks(1:n); xticklabels(ticklbls);
    xlim([0.9, n+0.5])
    title('CH validity index')
    legend({'Spatial','Spectral'})
    %%%
    ax3 = subplot(3,1,3); hold on
    yyaxis(ax3,'left')
    plot(kl_spat(1:n,1),'b*--');
    yyaxis(ax3,'right'), ax3.YColor = 'r';
    plot(kl_spec(1:n,1),'rx--');
    xticks(1:n); xticklabels(ticklbls);
    xlim([0.9, n+0.5])
    title('KL validity index')
    legend({'Spatial','Spectral'})
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % % for one column at a time % % %
crit_idx = 1; %5
sprintf('Dice = %0.3f \n DB = %0.2f / %0.2f \n DB_c = %0.2f / %0.2f \n time = %0.2f s', ...
    final_resuls_table(1,crit_idx), final_resuls_table(2,crit_idx), final_resuls_table(3,crit_idx), final_resuls_table(4,crit_idx), final_resuls_table(5,crit_idx), final_resuls_table(6,crit_idx))


% % % mean overall results % % %

mean_dice = nanmean(final_resuls_table(1,:))
mean_DBspat = nanmean(final_resuls_table(2,:))
mean_DBspec = nanmean(final_resuls_table(3,:))
mean_DBc_spat = nanmean(final_resuls_table(4,:))
mean_DBc_spec = nanmean(final_resuls_table(5,:))
mean_time = nanmean(final_resuls_table(6,:))
% nb_K

mean_res = nanmean(final_resuls_table,2)

