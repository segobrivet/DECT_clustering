


results_folder_name = "results_kmeans";
% results_folder_name = "results_sft";
% results_folder_name = "results_L75_K40";
% results_folder_name = "results_large_run";
pat_results = {};
nb_K = [];
dice_results = [];
% dice_mrg_flt_results = [];
time_results = [];
loglik_results = [];

% for patient_name = ["HNSCC2"	"HNSCC3"	"HNSCC5"	"HNSCC8"	"HNSCC9"	"HNSCC10"	"HNSCC11"	"HNSCC12"	"HNSCC13"	"HNSCC15"	"HNSCC15A"	"HNSCC17"	"HNSCC17A"	"HNSCC18"	"HNSCC20"	"HNSCC21"	"HNSCC22"	"HNSCC22A"	"HNSCC25"	"HNSCC26"]

patient_names = ["HNSCC2","HNSCC3","HNSCC5","HNSCC8","HNSCC9","HNSCC10",...
        "HNSCC11","HNSCC12","HNSCC13","HNSCC15","HNSCC15A","HNSCC17","HNSCC17A","HNSCC18","HNSCC20",...
        "HNSCC21","HNSCC22","HNSCC22A","HNSCC25","HNSCC26","HNSCC27","HNSCC29","HNSCC30",...
        "HNSCC31A","HNSCC32","HNSCC33","HNSCC34","HNSCC35","HNSCC36","HNSCC37A","HNSCC38","HNSCC39",...
        "HNSCC41","HNSCC42","HNSCC44","HNSCC44AM","HNSCC45","HNSCC46","HNSCC47","HNSCC48","HNSCC49",...
        "HNSCC51","HNSCC52","HNSCC52AM","HNSCC53","HNSCC55","HNSCC56","HNSCC57",...
        "HNSCC61A","HNSCC62","HNSCC63","HNSCC63A","HNSCC64A","HNSCC65A","HNSCC66","HNSCC67","HNSCC68","HNSCC69",...
        "HNSCC70A","HNSCC71","HNSCC72A","HNSCC73","HNSCC74","HNSCC75","HNSCC76","HNSCC77","HNSCC78","HNSCC79","HNSCC80",...
        "HNSCC81","HNSCC82","HNSCC83","HNSCC84","HNSCC85","HNSCC87","HNSCC88","HNSCC89","HNSCC90",...
        "HNSCC91","HNSCC92","HNSCC95","HNSCC96","HNSCC97","HNSCC98",...
        "HNSCC100","HNSCC101","HNSCC103","HNSCC105","HNSCC106","HNSCC108","HNSCC109"];
% 		% "HNSCC1","HNSCC10A","HNSCC60","HNSCC102"]


for patient_name = patient_names
    
%   lambda = [0.055,0.06,0.065,0.07,0.075,0.08,0.09,0.2]
    lambda = 0.075;
    roi_radius = 150;

%     fls = dir(fullfile(results_folder_name,[char(patient_name),'_*ROI',num2str(roi_radius),'_*',num2str(lambda),'_mixstats_red.mat']));
    fls = dir(fullfile(results_folder_name,[char(patient_name),'_*_imsegkmeans_mixstats_red.mat']));
    if length(fls) > 1
        disp("Several results are in folder for this case");
    elseif length(fls) == 0
        nb_K = [nb_K, NaN];
        dice_results = [dice_results, NaN];
        % dice_mrg_flt_results = [dice_mrg_flt_results, NaN];
        time_results = [time_results, NaN];
        % loglik_results = [loglik_results, NaN];
    else
        mixstatspath = fullfile(results_folder_name,fls(1).name);
        load(mixstatspath);
        
        nb_K = [nb_K, length(unique(mixstats_red.klas))];
        dice_results = [dice_results, mixstats_red.dice];
        % dice_mrg_flt_results = [dice_mrg_flt_results, mixstats_red.dice_mrg_flt];
        time_results = [time_results, mixstats_red.elaps_time];
        % loglik_results = [loglik_results, mixstats_red.loglik];
    end

end

mean_nbK = nanmean(nb_K)
std_nbK = nanstd(nb_K)
mean_dice = nanmean(dice_results)
std_dice = nanstd(dice_results)
mean_time = nanmean(time_results)
std_time = nanstd(time_results)
