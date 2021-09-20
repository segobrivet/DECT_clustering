
%% %%%%%%%% MIXMODEL TO SIGMA_K %%%%%%%% %%

% results_folder_name = 'results_3d_bigROIs_Fun_R07';
results_folder_name = 'results_old\results_3d_bigROIs_old\HNSCC9';

fls = dir(fullfile(results_folder_name,'/HNSCC9*cl50_mixmodel.mat'));
for nm = 1:length(fls)
    mixmodelpath = fullfile(results_folder_name,fls(nm).name);
    str = split(fls(nm).name,'_');
    patient_name = str{1}; rdn = str{2};
    roi_radius = str{4}; roi_radius = str2num(roi_radius(4:end));
    K = str{5}; K = str2num(K(3:end));
    fn_save_pdf = fullfile(results_folder_name, char(patient_name));

    load(mixmodelpath)
    mixmodel_sigma = mixmodel.Sigmak2;
    
    save([fn_save_pdf,'_',rdn,'__ROI',num2str(roi_radius),'_cl',num2str(K),'_sigmak2.mat'],'mixmodel_sigma')
    
end