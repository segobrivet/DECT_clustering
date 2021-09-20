
%% %%%%%%%% MIXSTATS TO LOGLIK %%%%%%%% %%

% results_folder_name = 'results_organs_Fun_R08_K';
% 
% fls = dir(fullfile(results_folder_name,'/subject*mixstats.mat'));
% for nm = 1:length(fls)
%     mixstatspath = fullfile(results_folder_name,fls(nm).name);
%     % mixmodelpath = [mixstatspath(1:end-12),'mixmodel.mat'];  % (1:end-17)
%     str = split(fls(nm).name,'_');
%     patient_name = str{1}; rdn = str{2};
%     roi_radius = str{4}; roi_radius = str2num(roi_radius(4:end));
%     K = str{5}; K = str2num(K(3:end));
%     fn_save_pdf = fullfile(results_folder_name, char(patient_name));
% 
%     load(mixstatspath)
%     loglik.curve = mixstats.stored_loglik;
%     
%     save([fn_save_pdf,'_',rdn,'__ROI',num2str(roi_radius),'_cl',num2str(K),'_loglik.mat'],'loglik')
%     
% end



%% %%%%%%%% PLOT LOGLIK %%%%%%%% %%

clear

pn = [];
gr = [];
negsum = [];
sumgr = [];
    
results_folder_name = 'results_organs_Fun_R07_K';

fls = dir(fullfile(results_folder_name,'/subject*cl40_loglik.mat'));
for nm = 1:length(fls)
    mixstatspath = fullfile(results_folder_name,fls(nm).name);
    % mixmodelpath = [mixstatspath(1:end-12),'mixmodel.mat'];  % (1:end-17)
    str = split(fls(nm).name,'_');
    patient_name = str{1}; rdn = str{2};
    roi_radius = str{4}; roi_radius = str2num(roi_radius(4:end));
    K = str{5}; K = str2num(K(3:end));
    fn_save_pdf = fullfile(results_folder_name, char(patient_name));

    close all
    
    load(mixstatspath)
    
    loglik.curve = stored_loglik;
    
    end_iter = length(loglik.curve);
    grd = gradient(loglik.curve(50:end_iter));
    sum_grd = 0;
    for i=1:length(grd)
        if grd(i) < 0
            sum_grd = sum_grd + grd(i);
        end
    end
    
    fig_loglik = figure('units','pixels','outerposition',[856   314   912   594]);
    plot(1:end_iter,loglik.curve,'b')
    hold on 
    plot(50:end_iter,grd,'r')
    legend(["loglik","gradient"],'Location','southeast')
    title(sprintf('Log-likelihood curve, end at l = %0.1f, sum of grad = %0.3f, sum of neg grad = %0.3f', loglik.curve(end_iter), sum(grd), sum_grd))
    
    save_pdf(fig_loglik,[fn_save_pdf,'_',rdn,'__ROI',num2str(roi_radius),'_cl',num2str(K),'_loglik']);
    
    loglik.grad = grd;
    loglik.negSumGrad = sum_grd;
    save([fn_save_pdf,'_',rdn,'__ROI',num2str(roi_radius),'_cl',num2str(K),'_loglik'],'loglik');
     
    pn = [pn, string(patient_name)];
    gr = [gr, loglik.curve(end)];
    negsum = [negsum, sum_grd];
    sumgr = [sumgr, sum(grd)];
    
end


pn
gr
negsum
sumgr
