

results_folder_name = 'results_3d_bigROIs_R07';
% results_folder_name = 'results_3d_bigROIs_Rdiag033';
% results_folder_name = 'results_test';

fls = dir(fullfile(results_folder_name,'/HNS*stats.mat'));

for nm = 1:length(fls)
    
    load(fullfile(results_folder_name,fls(nm).name));
    
    str = split(fls(nm).name,'_');
    patient_name = str{1};
    rdn = str{2};
    roi_radius = str{4};
    roi_radius = str2num(roi_radius(4:end));
    K = str{5};
    K = str2num(K(3:end));
    
    mixstats_klas = mixstats.klas;
    
    fn_save_pdf = fullfile(results_folder_name, char(patient_name));
    save([fn_save_pdf,'_',rdn,'__ROI',num2str(roi_radius),'_cl',num2str(K),'_mixstats_klas.mat'],'mixstats_klas')
    
end