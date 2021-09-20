%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Classifying curves
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Choose specifications

% clear all; close all;
% clc;

addpath('utils');

%%%%%%%%%%%%%%%%%%%%%%%%%%% MACHINE and FOLDER NAMES %%%%%%%%%%%%%%%%%%%%%%
% machine_type = 'sego_laptop';
machine_type = 'sego_clean';
% machine_type = 'Peter';
% machine_type = 'GPU';

traindata_folder_name = 'results_clustering_L07_K60';

testdata_folder_name = 'results_clustering_L07_K60';

results_folder_name = 'results_classif_L07_K60';
% results_folder_name = 'ML_results_organ_decayCurves';


%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% classify curves or clusters
type_of_data = 'curves';
% type_of_data = 'clusters';

% ROI size
roi_radius = 90;

% Nb of clusters
K = 60;

% tissue enhancement window
lvl = 150;  % 40  % 70
wdw = 700;  % 350  % 500

% organs ID
org_ids = '1-2';


%%%%%%%%%%%%%%%%%%%%%%%%% Model Fit or Load %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% neuralnetwork = 'nn';
neuralnetwork = 'lstm';
% neuralnetwork = 'ML';

run_train = 1;
load_saved_mdl = 0;

run_test = 1;
do_initResults = 1;
load_cpred = 0;
do_clusterResults = 1;
load_propTrue = 0;

save_training = 1;
save_results = 1;
rdm = num2str(randi(1000));


%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTS OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
show_slices = 'middle';
% show_slices = 'bottom';
% show_slices = 'top';
max_slic_subplot = 6;

plot_trainSubjects = 1;
plot_results = 1;



if run_train
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                            TRAINING DATA                                %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    disp('Train run');
    
    % Build dataset
    YTr = [];
    CTr = [];
    pat_concat = [];
    if plot_trainSubjects, fig_curves = figure('units','normalized','outerposition',[0 0 1 1]); spidx = 0; end
    
    %%%%%%%%%%
    % Run on multiple patients from a folder
    % fls = dir(fullfile(traindata_folder_name,'/HNSCC*mixstats.mat'));
%     fls = dir(fullfile(traindata_folder_name,'/subject_*mixstats_red.mat'));
%     % fls = dir(fullfile(traindata_folder_name,'/HNSCC*cl50_mixstats_red.mat'));
%     nb_patients = length(fls)-1;
%     for nm = 1:nb_patients
%         str = split(fls(nm).name,'_');
%         patient_name = str{1}; rdn = str{2};
%         roi_radius = str{4}; roi_radius = str2num(roi_radius(4:end));
%         K = str{5}; K = str2num(K(3:end));
%         mixstatspath = fullfile(traindata_folder_name,fls(nm).name);
%         mixmodelpath = [mixstatspath(1:end-12),'mixmodel.mat']; % from mixstats

%     fl_vec = [%"results_3d_bigROIs_Fun_R07/HNSCC9_585__ROI90_cl50_mixstats_red.mat",
%               "results_3d_bigROIs_extraFeat_R09_cl70/HNSCC9_126__ROI90_cl70_mixstats_red.mat",
%               "results_3d_bigROIs_extraFeat_R10/HNSCC9_65__ROI90_cl50_mixstats_red.mat"];
%     fl_vec = ["results_test/HNSCC9_944__ROI20_cl7_mixstats.mat",
%               "results_3d_bigROIs_Fun_R07/HNSCC9_585__ROI90_cl50_mixstats_red.mat"];
%     nb_patients = length(fl_vec);
%     for nm = 1:nb_patients
%         fl_nm = split(fl_vec(nm),'/');
%         str = split(fl_nm(end),'_');
%         patient_name = str{1}; rdn = str{2};
%         roi_radius = str{4}; roi_radius = str2num(roi_radius(4:end));
%         K = str{5}; K = str2num(K(3:end));
%         mixstatspath = char(fullfile(join(fl_nm(1:2),'/')));
%         mixmodelpath = [mixstatspath(1:end-12),'mixmodel.mat']; % from mixstats
% %         sigmak2path = [mixstatspath(1:end-16),'sigmak2.mat'];  % from mixstats_red
    
    %%%%%%%%%% OR %%%%%%%%%%
    
    % % Run on multiple patients manually
%     patient_names = ["subject3","subject18","subject21","subject23","subject8","subject15","subject94"];
    % patient_names = ["subject3","subject12","subject14","subject18","subject21","subject23","subject8","subject15","subject31","subject93","subject94"];

%         PB? "HNSCC1",
%         TO BUILD: "HNSCC20","HNSCC10A"
%           check slices<->GT on subj 15A, 17A, 61A?
% 	patient_names = ["HNSCC22","HNSCC15A","HNSCC17A","HNSCC22A","HNSCC44AM","HNSCC52AM","HNSCC61A","HNSCC63A","HNSCC64A","HNSCC65A","HNSCC27","HNSCC29","HNSCC30","HNSCC31A","HNSCC33","HNSCC35","HNSCC36","HNSCC37A","HNSCC38","HNSCC39","HNSCC41","HNSCC42","HNSCC44","HNSCC45","HNSCC46","HNSCC47","HNSCC51","HNSCC52","HNSCC53","HNSCC55","HNSCC56","HNSCC57","HNSCC60","HNSCC66","HNSCC67","HNSCC69"];
    % 1st round % ["subject8_tumor","HNSCC2","HNSCC3","HNSCC5","HNSCC8","HNSCC9","HNSCC10","HNSCC11","HNSCC12","HNSCC13","HNSCC15","HNSCC17","HNSCC18","HNSCC26"];
    % 2nd round % ["HNSCC1","HNSCC20","HNSCC22","HNSCC27","HNSCC29","HNSCC30","HNSCC31A","HNSCC33","HNSCC35","HNSCC36","HNSCC37A","HNSCC38","HNSCC39","HNSCC41",...
                % "HNSCC42","HNSCC44","HNSCC45","HNSCC46","HNSCC47","HNSCC51","HNSCC52","HNSCC53","HNSCC55","HNSCC56","HNSCC57","HNSCC60","HNSCC66","HNSCC67","HNSCC69",...
                % "HNSCC10A","HNSCC15A","HNSCC17A","HNSCC22A","HNSCC44AM","HNSCC52AM","HNSCC61A","HNSCC63A","HNSCC64A","HNSCC65A"]

    patient_names = ["HNSCC63A"];
    nb_patients = length(patient_names);
    for patient_name = patient_names
       
        %%%%%%%%%%
        
        patient_nm = char(patient_name); disp(patient_nm);
        pat_concat = [patient_nm,'_',pat_concat];
        
        % Build dataset for a specific patient
        take_all_slices = 0;
        res = build_patient_data(machine_type, patient_nm, take_all_slices, org_ids, roi_radius, max_slic_subplot, show_slices);
        
        Y = res.Y;
        decay_curves = res.decay_curves;
        klas_tum = res.klas_tum;
        lin_obj = res.lin_obj; gr_truth = res.gr_truth;
        rmin = res.rmin; rmax = res.rmax; cmin = res.cmin; cmax = res.cmax;
        slic_show = res.slic_show; subj_slic = res.subj_slic; slic_min_idx = res.slic_min_idx;
        % take a margin around boundary
        klas_tum = ismember(lin_obj, find(imerode(gr_truth,strel('disk',3))));
        klas_obj = ~ismember(lin_obj, find(imdilate(gr_truth,strel('disk',3))));
        
        switch(type_of_data)
            
            case('curves')
                % Balance dataset proportions between true (organs, tumor...) and false (rest of body)
                ind_obj = find(klas_obj);
                ind_rand = randperm(length(ind_obj),length(ind_obj)-sum(klas_tum));
                ind_del = [ind_obj(ind_rand); find(~(klas_tum + klas_obj))];
                Y(ind_del,:) = [];
                decay_curves(ind_del,:) = [];
                klas_tum(ind_del) = [];
                lin_obj(ind_del) = [];

                % Concatenate multiple patients datasets
%                 YTr = [YTr; Y];
                YTr = [YTr; decay_curves];
                CTr = [CTr; klas_tum];
                
                if plot_trainSubjects
                    % Plot training voxels of current patient
                    reco_lbl = -ones(size(gr_truth));
                    reco_lbl(lin_obj) = klas_tum;
                    for i=1:length(slic_show)
                        tumor_contour_list{i} = bwboundaries(gr_truth(:,:,slic_show(i)));   % shift of -1 after slic_show(i)?
                    end
                    fig_slic = plot_classif_results(reco_lbl, tumor_contour_list, subj_slic, slic_show, max_slic_subplot, slic_min_idx, lvl, wdw, rmin, rmax, cmin, cmax);
                    suptitle(sprintf('6 %s slices / %d.   Voxels involved in training (on last loaded patient)', show_slices, size(subj_slic,3)));

%                     figure(fig_curves);
%                     spidx = spidx + 1;
%                     subplot(2,nb_patients,spidx)
%                     plot(1:21, decay_curves(find(klas_tum==0),:))
%                     xlabel('keV'), ylabel('HU')
%                     xlim([1 21]), ylim([-1000 3200])
%                     title([patient_nm,' - ',num2str(sum(klas_tum==0)),' body curves'])
%                     subplot(2,nb_patients,nb_patients+spidx)
%                     plot(1:21, decay_curves(find(klas_tum),:))
%                     xlabel('keV'), ylabel('HU')
%                     xlim([1 21]), ylim([-500 1000])
%                     title([patient_nm,' - ',num2str(sum(klas_tum)),' tumor curves'])
%                     if spidx == nb_patients, sgtitle("Energy decay curves used in training",'Fontweight','bold'); end
                    
                end
                

            case('clusters')
                load(mixstatspath)
                load(mixmodelpath)
%                 load(mixstatspath)
%                 load(sigmak2path)
%                 mixstats = mixstats_red;
%                 mixmodel.Sigmak2 = mixmodel_sigma;
                
                Y = cell(K,1);
%                 for k=1:K
%                     Y{k} = [mixstats.Muk(:,k)';
%                             mixstats.Muk(:,k)'-2*sqrt(mixmodel.Sigmak2(k));
%                             mixstats.Muk(:,k)'+2*sqrt(mixmodel.Sigmak2(k))];
%                 end

                % Compute features
                % % TODO
                
                match_klas = match_tumor_clusters(klas_tum,mixstats.klas,K);
                C = zeros(K,1);
                C(match_klas) = 1;
                
                YTr = [YTr; Y];
                CTr = [CTr; C];
                
                
                if plot_trainSubjects
                if K < 64
                    % Plot normalized curves
                    fig_normCurves = figure('units','normalized','outerposition',[0 0 0.5 0.7]);
                    ha = set_subplot_grid(K);
                    for k=1:K %min(K,63)
                        sigmak2 = sqrt(mixmodel.Sigmak2(k));
                        axes(ha(k));
                        if ismember(k,match_klas)
                            plot(1:size(YTr{1},2),mixstats.Muk(:,k)','color',[1 0 0],'linewidth',5);
                            hold on
                            plot(1:size(YTr{1},2),mixstats.Muk(:,k)-2*sigmak2,'--r','linewidth',1)
                            plot(1:size(YTr{1},2),mixstats.Muk(:,k)+2*sigmak2,'--r','linewidth',1)
                            ylabel(['TUMOR - Cl ',num2str(k)],'FontWeight','bold');
                        else
                            plot(1:size(YTr{1},2),mixstats.Muk(:,k)','color',[0 0 1],'linewidth',5);
                            hold on
                            plot(1:size(YTr{1},2),mixstats.Muk(:,k)-2*sigmak2,'--b','linewidth',1)
                            plot(1:size(YTr{1},2),mixstats.Muk(:,k)+2*sigmak2,'--b','linewidth',1)
                            ylabel(['Cluster ',num2str(k)],'FontWeight','bold');
                        end
                        % xlim([min(T) max(T)]);
                        % ylim([min(Y,[],'all') max(Y,[],'all')]);
                    end
                    box on;
                    suptitle(['Normalized curves per cluster'])
                    for k = (K+1):length(ha)
                        axes(ha(k));
                        axis off
                    end
                end
                end
    
    
            otherwise
                %
        end
        
    end

    
    % % % % % %
    
    fn_save_mdl = fullfile(results_folder_name, pat_concat);  % to save fig as pdf
    
    if plot_trainSubjects
%         figure;
%         plot(1:21, YTr(find(CTr==0),:))
%         title('All training curves from body')
%         figure;
%         plot(1:21, YTr(find(CTr),:))
%         title('All training curves from tumor')
        
        % save image in pdf
        if  save_training
%             savefig(fig_slic,[fn_save_mdl,'_',rdm,'_trainingVoxels']);
            save_pdf(fig_slic, [fn_save_mdl,'_',rdm,'_trainingVoxels.pdf']);
%             save_pdf(fig_curves, [fn_save_mdl,'_',rdm,'_trainingCurves.pdf']);
        end
        
    end
    
    % % % % % %
    
    
    % Split Train-Val dataset
    prop_val = 0.2;
    cvtv = cvpartition(size(YTr,1),'HoldOut',prop_val);
    YTrain = YTr(~cvtv.test,:);
    CTrain = CTr(~cvtv.test);
    YVal  = YTr(cvtv.test,:);
    CVal = CTr(cvtv.test);
    
    %     YTrain = Y;
    %     CTrain = klas_tum;
    %
    %     ind_obj = find(CTrain == 0);
    %     ind_rand = randperm(length(ind_obj),length(CTrain)-2*(length(CTrain)-length(ind_obj)));
    %     YTrain(ind_obj(ind_rand),:) = [];
    %     CTrain(ind_obj(ind_rand)) = [];
    %

    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                            MODEL FITTING                                %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Bi-LSTM FITTING
    
    if load_saved_mdl
        
        fls = dir(fullfile(results_folder_name,'*815_ML_model.mat'));   % tum: 815, org:
        model_path = fullfile(results_folder_name,fls(1).name);
        
        load(model_path)
        
    else
        disp('Starting learning...'); tic;
        
        switch(neuralnetwork)
            
            case('svm')
                
            case('nn')
                
                net = train(fitnet([10 3]),YTrain',CTrain');
                
                CPred = net(YTrain');
                perf = perform(net,CPred,CTrain')
                
            case('lstm')
                
                numFeatures = 1;
                numHiddenUnits = 100;
                numClasses = 2;
                
                layers = [ ...
                    sequenceInputLayer(numFeatures)
                    bilstmLayer(numHiddenUnits,'OutputMode','sequence') %
                    dropoutLayer(0.2) %
                    bilstmLayer(numHiddenUnits,'OutputMode','last')
                    fullyConnectedLayer(numClasses)
                    softmaxLayer
                    classificationLayer]
                
                maxEpochs = 8;
                miniBatchSize = 1024;
                
                options = trainingOptions('adam', ...  % sgdm
                    'ExecutionEnvironment','auto', ...
                    'MaxEpochs',maxEpochs, ...
                    'MiniBatchSize',miniBatchSize, ...
                    'LearnRateSchedule','piecewise', ...
                    'ValidationData',{num2cell(YVal,2),categorical(CVal)}, ...
                    'Shuffle','every-epoch', ... % once
                    'Verbose',1, ...
                    'Plots','training-progress');
%                     'GradientThreshold',1, ...
%                     'LearnRateDropFactor',0.2, ...
%                     'LearnRateDropPeriod',5, ...
                
                net = trainNetwork(num2cell(YTrain,2),categorical(CTrain),layers,options);
                
            otherwise
                %
        end
        fprintf('Elapsed time %f sec \n', toc);
    end
    
    % % % % % %
    
    % Trainning accuracy
    CPredTrain = classify(net,num2cell(YTrain,2), 'SequenceLength','longest');
    cpredtrain = string(CPredTrain) == '1';
    acc_train = sum(cpredtrain == CTrain)./numel(CTrain)
    
    CPredVal = classify(net,num2cell(YVal,2), 'MiniBatchSize',miniBatchSize, 'SequenceLength','longest');
    cpredval = string(CPredVal)=='1';
    acc_val = sum(cpredval == CVal)./numel(CVal)
        
    % % % % % %
    
    % Save model
    if save_training  &&  ~load_saved_mdl
        save([fn_save_mdl,'_',rdm,'_ML_model.mat'],'net')
    end
    
    % % % % % %
    
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               RESULTS                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Build test data from 1 patient
if run_test
    
    disp('Test run');
    subj_test = 'subject12'; disp(subj_test);  % H9(just all)  % on subj12, 92, 14, 31, 93, 24
    
    if load_saved_mdl
        fls = dir(fullfile(results_folder_name,'*815_ML_model.mat'));
        model_path = fullfile(results_folder_name,fls(1).name);
        load(model_path)
    end
    
    % Choose patient
    fls = dir(fullfile(testdata_folder_name,[subj_test,'*cl50_mixstats_red.mat']));
    nm = length(fls);
    mixstatspath = fullfile(testdata_folder_name,fls(nm).name);
    str = split(fls(nm).name,'_');
    patient_nm = char(str{1}); rdn = str{2};
    roi_radius = str{4}; roi_radius = str2num(roi_radius(4:end));
    K = str{5}; K = str2num(K(3:end));
    fn_save_results = fullfile(results_folder_name, patient_nm);  % to save fig as .fig and .pdf
    
    
    % Build dataset for this patient
    take_all_slices = 0; % 0 = only 8 middle tumor slices, 1 = all tumorous slices
    res = build_patient_data(machine_type, patient_nm, take_all_slices, org_ids, roi_radius, max_slic_subplot, show_slices);
    lin_obj = res.lin_obj;
    gr_truth = res.gr_truth;
    rmin = res.rmin; rmax = res.rmax; cmin = res.cmin; cmax = res.cmax;
    slic_show = res.slic_show; subj_slic = res.subj_slic; slic_min_idx = res.slic_min_idx;
    coord = res.coord;
    
%     YTest = res.Y;
    YTest = res.decay_curves;
    CTest = res.klas_tum;
    
    for i=1:length(slic_show)
        tumor_contour_list{i} = bwboundaries(gr_truth(:,:,slic_show(i)));   % shift of -1 after slic_show(i)?
    end
    
    
    
    if do_initResults
        
        if load_cpred
            
            flp = dir(fullfile(results_folder_name,[subj_test,'*_cpredtestall.mat']));   % 815
            cpred_path = fullfile(results_folder_name,flp(1).name);
            load(cpred_path)
            
        else
            
            % Classify
            CPred = classify(net,num2cell(YTest,2), 'SequenceLength','longest');
            cpred = string(CPred)=='1';
            acc_test = sum(cpred == CTest)./numel(CTest)
            
        end
        
        % % % %
        
        if plot_results
            
            % figure;
            % cm = confusionchart(categorical(CTest),categorical(cpred));
            figure;
            plotconfusion(categorical(CTest),categorical(cpred));

            % % % %

            % Reconstruct results as 3D volume
            reco_lbl = -ones(size(gr_truth));
            reco_lbl(lin_obj) = cpred;

            % Plot results
            fig_slic = plot_classif_results(reco_lbl, tumor_contour_list, subj_slic, slic_show, max_slic_subplot, slic_min_idx, lvl, wdw, rmin, rmax, cmin, cmax);
            suptitle(sprintf('6 %s slices / %d.   Classification test results', show_slices, size(subj_slic,3)));

            % % % %

            % save image in pdf
            if save_results
                savefig(fig_slic,[fn_save_results,'_',rdm,'_classifResults_all']);
                save_pdf(fig_slic, [fn_save_results,'_',rdm,'__classifResults_all.pdf']);
            end
        
        else
            
            disp('Confusion matrix:');
            disp(confusionmat(categorical(CTest),categorical(cpred)));
            
            if save_results
                save([fn_save_results,'_',rdm,'_cpredtestall.mat'],'cpred')
            end
        end
        
        
    end
    
    
    if do_clusterResults
        
        % load clustering results
        load(mixstatspath);
        mixstats.klas = mixstats_red.klas;
        
        % prop_true records the proportion of tumor voxel among all voxels in each cluster
        if load_propTrue
            flp = dir(fullfile(results_folder_name,[subj_test,'*_propTrue_K.mat']));   % 815
            propTrue_path = fullfile(results_folder_name,flp(1).name);
            load(propTrue_path)
        else
            prop_true = -ones(1,K);
            for cl_id=1:K
                ind_cl = mixstats.klas == cl_id;
                switch(neuralnetwork)
                    case('nn')
                        CPred = net(YTest');
                        mean_cpred = mean(CPred == 1);
                    case('lstm')
                        CPred_k = classify(net,num2cell(YTest(ind_cl,:),2), 'SequenceLength','longest');
                        mean_cpred = mean(string(CPred_k) == '1'); % sum(cpred)
                    otherwise
                        %
                end
                prop_true(cl_id) = mean_cpred;
            end
        end
        
        
        if plot_results
            
            % % % Determine threshold value to select tumor clusters
            % Select the cluster with the highest proportion
            [sorted_prop, sorted_idx] = sort(prop_true,'descend');
            si = 1;
            thresh = sorted_prop(si);
            select_coord = coord(mixstats.klas == sorted_idx(si),:);
            % aggregate clusters one-by-one, only if it is neighboring the selection, stop otherwise
            for si = 2:K
                % point cloud coordinates of considered cluster
                curr_coord = coord(mixstats.klas == sorted_idx(si),:);
                % closest point of considered-cluster-point-cloud to the mean of already-selected-cluster-point-cloud
                closest_pt_idx = knnsearch(curr_coord,mean(select_coord,1));
                closest_point = curr_coord(closest_pt_idx,:);
                % closest point of already-selected-cluster-point-cloud to this previous point-cloud-border point
                closest_sl_idx = knnsearch(select_coord,closest_point);
                closest_select = select_coord(closest_sl_idx,:);
                % compute cluster distance
                cluster_dist = sum((closest_point - closest_select).^2).^0.5;
                
                if cluster_dist < 8   % 8
                    thresh = sorted_prop(si);
                    select_coord = [select_coord; curr_coord];
                else
                    break;
                end
            end
            % create result image
            reco_lbl = -ones(size(gr_truth));
            for ii = 1:(si-1)
                ind_cl = mixstats.klas == sorted_idx(ii);
                reco_lbl(lin_obj(ind_cl)) = 1;
            end
            for ii=si:K
                ind_cl = mixstats.klas == sorted_idx(ii);
                reco_lbl(lin_obj(ind_cl)) = 0;
            end
            
            %%% OLD %%%
%             thresh = 0.75;
%             % create result image
%             reco_lbl = -ones(size(gr_truth));
%             for cl_id=1:K
%                 ind_cl = mixstats.klas == cl_id;
%                 reco_lbl(lin_obj(ind_cl)) = prop_true(cl_id) >= thresh;
%             end
            %%% END OLD %%%
            
            % Plot classification results
            fig_slic2 = plot_classif_results(reco_lbl, tumor_contour_list, subj_slic, slic_show, max_slic_subplot, slic_min_idx, lvl, wdw, rmin, rmax, cmin, cmax);
            suptitle(sprintf('6 %s slices / %d.   Classification test results', show_slices, size(subj_slic,3)));

            % Cluster selection from tumor proportions
            fig_hist = figure();
            bar(prop_true)
            hold on
            plot([0 K+1],[thresh thresh], '--', 'LineWidth', 2, 'color', [0.8 0.3 0.1])
            title('Proportion of tumor voxels per cluster')
            xlabel('Cluster Id')
            ylabel('Proportion')
            yticks(0:0.05:1)
            grid on

            % save image in pdf
            if save_results
                savefig(fig_slic2,[fn_save_results,'_',rdm,'_classifResults_K']);
                save_pdf(fig_slic2, [fn_save_results,'_',rdm,'__classifResults_K.pdf']);
                save([fn_save_results,'_',rdm,'__propTrue_K.mat'],'prop_true');
                save_pdf(fig_hist, [fn_save_results,'_',rdm,'__hist_propTrue.pdf']);
            end
            
        else
            
            if save_results
                save([fn_save_results,'_',rdm,'__propTrue_K.mat'],'prop_true');
            end
            
        end
        
        
    end
    
    
end




