
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                           Analyse all results                           %
%                                                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To run this script, you first need to:
% 1. main.m: learn model and get clustering results in 'mixstats_red' variable
% 2. build_cluster_idx.m: compute clustering index scores and Dice scores
% 3. build_results_table.m: crun and save results as .mat file
% 4. once results are saved for every method, come here and launch analysis


load('results_SgMFR-Bspl.mat')
results_SgMFR_Bspl = results;
load('results_SgMFR-poly.mat')
results_SgMFR_poly = results;
load('results_SgMVFR-Bspl.mat')
results_SgMVFR_Bspl = results;
load('results_SsMFR-Bspl.mat')
results_SsMFR_Bspl = results;
load('results_GMM.mat')
results_GMM = results;
load('results_kmeans.mat')
results_kmeans = results;
load('results_SelSearch.mat')
results_SelSearch = results;


%% T-test
% Significance of improvement between best of baseline and worse of our method.
[h,p] = ttest2(results_SgMFR_poly.dice,results_kmeans.dice)



%% BOXPLOT for metrics with all results

figure; hold on

% DB spatial
subplot(1,5,1)
hold on
scatter(1*ones(size(results_GMM.DB_spat)),results_GMM.DB_spat,'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
scatter(2*ones(size(results_kmeans.DB_spat)),results_kmeans.DB_spat,'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
scatter(3*ones(size(results_SelSearch.DB_spat)),results_SelSearch.DB_spat,'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
scatter(4*ones(size(results_SgMFR_Bspl.DB_spat)),results_SgMFR_Bspl.DB_spat,'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
scatter(5*ones(size(results_SgMFR_poly.DB_spat)),results_SgMFR_poly.DB_spat,'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
scatter(6*ones(size(results_SgMVFR_Bspl.DB_spat)),results_SgMVFR_Bspl.DB_spat,'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
scatter(7*ones(size(results_SsMFR_Bspl.DB_spat)),results_SsMFR_Bspl.DB_spat,'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);

boxplot([results_GMM.DB_spat',results_kmeans.DB_spat',results_SelSearch.DB_spat',...
    results_SgMFR_Bspl.DB_spat',results_SgMFR_poly.DB_spat',results_SgMVFR_Bspl.DB_spat',results_SsMFR_Bspl.DB_spat'])

xticklabels(["GMM","k-means","S.Search","SgMFR-Bspl","SgMFR-poly","SgMVFR-Bspl","SsMFR-Bspl"]), xtickangle(90)
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
title('DB - Spatial index')

% DB spectral
results_SsMFR_Bspl.DB_spec(results_SsMFR_Bspl.DB_spec > 30) = 25;
subplot(1,5,2)
hold on
scatter(1*ones(size(results_GMM.DB_spec)),results_GMM.DB_spec,'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
scatter(2*ones(size(results_kmeans.DB_spec)),results_kmeans.DB_spec,'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
scatter(3*ones(size(results_SelSearch.DB_spec)),results_SelSearch.DB_spec,'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
scatter(4*ones(size(results_SgMFR_Bspl.DB_spec)),results_SgMFR_Bspl.DB_spec,'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
scatter(5*ones(size(results_SgMFR_poly.DB_spec)),results_SgMFR_poly.DB_spec,'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
scatter(6*ones(size(results_SgMVFR_Bspl.DB_spec)),results_SgMVFR_Bspl.DB_spec,'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
scatter(7*ones(size(results_SsMFR_Bspl.DB_spec)),results_SsMFR_Bspl.DB_spec,'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);

boxplot([results_GMM.DB_spec',results_kmeans.DB_spec',results_SelSearch.DB_spec',...
    results_SgMFR_Bspl.DB_spec',results_SgMFR_poly.DB_spec',results_SgMVFR_Bspl.DB_spec',results_SsMFR_Bspl.DB_spec'])

xticklabels(["GMM","k-means","S.Search","SgMFR-Bspl","SgMFR-poly","SgMVFR-Bspl","SsMFR-Bspl"]), xtickangle(90)
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
title('DB - Spectral index')

% DB for tumor spatial
results_SsMFR_Bspl.DBc_spat(results_SsMFR_Bspl.DBc_spat > 80) = 75;
subplot(1,5,3)
hold on
scatter(1*ones(size(results_GMM.DBc_spat)),results_GMM.DBc_spat,'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
scatter(2*ones(size(results_kmeans.DBc_spat)),results_kmeans.DBc_spat,'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
scatter(3*ones(size(results_SelSearch.DBc_spat)),results_SelSearch.DBc_spat,'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
scatter(4*ones(size(results_SgMFR_Bspl.DBc_spat)),results_SgMFR_Bspl.DBc_spat,'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
scatter(5*ones(size(results_SgMFR_poly.DBc_spat)),results_SgMFR_poly.DBc_spat,'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
scatter(6*ones(size(results_SgMVFR_Bspl.DBc_spat)),results_SgMVFR_Bspl.DBc_spat,'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
scatter(7*ones(size(results_SsMFR_Bspl.DBc_spat)),results_SsMFR_Bspl.DBc_spat,'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);

boxplot([results_GMM.DBc_spat',results_kmeans.DBc_spat',results_SelSearch.DBc_spat',...
    results_SgMFR_Bspl.DBc_spat',results_SgMFR_poly.DBc_spat',results_SgMVFR_Bspl.DBc_spat',results_SsMFR_Bspl.DBc_spat'])

xticklabels(["GMM","k-means","S.Search","SgMFR-Bspl","SgMFR-poly","SgMVFR-Bspl","SsMFR-Bspl"]), xtickangle(90)
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
title('tumor DB - Spatial index')

% DB for tumor spectral
subplot(1,5,4)
hold on
scatter(1*ones(size(results_GMM.DBc_spec)),results_GMM.DBc_spec,'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
scatter(2*ones(size(results_kmeans.DBc_spec)),results_kmeans.DBc_spec,'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
scatter(3*ones(size(results_SelSearch.DBc_spec)),results_SelSearch.DBc_spec,'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
scatter(4*ones(size(results_SgMFR_Bspl.DBc_spec)),results_SgMFR_Bspl.DBc_spec,'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
scatter(5*ones(size(results_SgMFR_poly.DBc_spec)),results_SgMFR_poly.DBc_spec,'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
scatter(6*ones(size(results_SgMVFR_Bspl.DBc_spec)),results_SgMVFR_Bspl.DBc_spec,'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
scatter(7*ones(size(results_SsMFR_Bspl.DBc_spec)),results_SsMFR_Bspl.DBc_spec,'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);

boxplot([results_GMM.DBc_spec',results_kmeans.DBc_spec',results_SelSearch.DBc_spec',...
    results_SgMFR_Bspl.DBc_spec',results_SgMFR_poly.DBc_spec',results_SgMVFR_Bspl.DBc_spec',results_SsMFR_Bspl.DBc_spec'])

xticklabels(["GMM","k-means","S.Search","SgMFR-Bspl","SgMFR-poly","SgMVFR-Bspl","SsMFR-Bspl"]), xtickangle(90)
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
title('tumor DB - Spectral index')


% Dice
subplot(1,5,5)
hold on
scatter(1*ones(size(results_GMM.dice)),results_GMM.dice,'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
scatter(2*ones(size(results_kmeans.dice)),results_kmeans.dice,'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
scatter(3*ones(size(results_SelSearch.dice)),results_SelSearch.dice,'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
scatter(4*ones(size(results_SgMFR_Bspl.dice)),results_SgMFR_Bspl.dice,'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
scatter(5*ones(size(results_SgMFR_poly.dice)),results_SgMFR_poly.dice,'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
scatter(6*ones(size(results_SgMVFR_Bspl.dice)),results_SgMVFR_Bspl.dice,'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
scatter(7*ones(size(results_SsMFR_Bspl.dice)),results_SsMFR_Bspl.dice,'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);

boxplot([results_GMM.dice',results_kmeans.dice',results_SelSearch.dice',...
    results_SgMFR_Bspl.dice',results_SgMFR_poly.dice',results_SgMVFR_Bspl.dice',results_SsMFR_Bspl.dice'])

xticklabels(["GMM","k-means","S.Search","SgMFR-Bspl","SgMFR-poly","SgMVFR-Bspl","SsMFR-Bspl"]), xtickangle(90)
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
title('Dice similarity score')

% save_pdf(gcf,['boxplots_results_all.pdf']);
