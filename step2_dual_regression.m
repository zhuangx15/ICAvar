function step2_dual_regression(freq_f1,loc_f1_Yf,freq_f1_date,ica_save_path,iq_threshold,...
    zscore_thre,interval_days)
% Input:
%   freq_f1: NsamplexNlocation; iVar produced alternative allele frequencies at
%   eahc mutation location;
%   freq_f1_date: Nsample x 1; contains sampling date;
%   loc_f1_Yf: Nlocation x 1; mutation location in freq_f1
%   ica_save_path: contains icasso/ica results;
%   iq_threshold: Estimate stability index; [0,1], 1--> most stable;
%   zscore_thre: zscore-threshold on ICA source file; 2 --> keep 5 positions (following Gaussian distribution);
%   interval_days: interval days to group wastewater samples together for strain detection through annotate dual-regressed source matrix
% Output:
%   dual-regressed source file for each interval time period stored under [ica_save_path,'\dual_regression'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist("interval_days",'var')
    interval_days = 7;
end
if ~exist("zscore_thre",'var')
    zscore_thre = 2;
end
if ~exist("iq_threshold",'var')
    iq_threshold = 0.7;
end
if ~exist('freq_f1','var') || ~exist("freq_f1_date",'var')
    error('Please provide sample files');
end
if ~exist('ica_save_path','var')
    error('Please specify ica run path');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ica_results_file = [ica_save_path,'\ica_results.mat'];
if ~exist(ica_results_file,'file')
    error('No ICA results file found, please run step 1 first');
end

freq_f1_full = freq_f1;
load([ica_save_path,'\ica_results.mat'],'sR','S','A','loc_f1','iq');
ind_ica_keep = find(iq>=iq_threshold);  % the estimates are concentrated in m compact and close-to-orthogonal clusters.
                                        % In this case the index to all estimate-clusters is (very close) to one.
                                        % The value drops when the clusters grow wider and mix up
                                        % selected by looking at similarity graph h4
loc_f1_ica = loc_f1;
%%%=================match positions in ICA and in data==========
[loc_f1,ind11,ind22] = intersect(loc_f1_ica,loc_f1_Yf,'stable');
Nloc_f1 = numel(loc_f1);
S_full = S(:,ind11);
freq_f1_full = freq_f1_full(:,ind22);
%%================================================================
S = S_full(ind_ica_keep,:);
nica = size(S,1);
S_zscore = zeros(size(S));
for j = 1:nica
    S_zscore(j,:) = zscore(S(j,:));
end
S_threshold = abs(S_zscore) >= zscore_thre; %threshold_for_significant_SNPs_on_S;
h7 = figure(18852);
h7 = signalplot(S_threshold,h7);
xlabel('Genome-position');
ylabel('Source signal (zscore)');
title(['Thresholded at ',num2str(zscore_thre)]);
set(gca,'FontSize',16);
box on; grid on;
saveas(h7,[ica_save_path,'\S_threshold@',num2str(zscore_thre),'.jpg']);

%% dual regression [might not be right] since A is already the loading for that sample
monthly_type = ['dual_regression_by_',num2str(interval_days),'days'];  %dual_regression_by_month
freq_f1_sample_date_num = datenum(freq_f1_date,'mmddyy');
min_date_num = min(freq_f1_sample_date_num);
max_date_num = max(freq_f1_sample_date_num);

date_vec_num = min_date_num:interval_days:max_date_num;
date_vec = cellstr(datestr(date_vec_num,'mm-dd-yyyy'));
dual_regression_save_path = [ica_save_path,'\',monthly_type];
if ~exist(dual_regression_save_path,'dir')
    mkdir(dual_regression_save_path);
end
Ndate = numel(date_vec);

S_sub = nan(Ndate,nica,Nloc_f1);
S_sub_zscore = nan(Ndate,nica,Nloc_f1);
Nsamplet = nan(Ndate,1);
for id_date = 1:Ndate
    ind_sample_tmp = find((freq_f1_sample_date_num-date_vec_num(id_date))<interval_days & (freq_f1_sample_date_num-date_vec_num(id_date))>=0);
    if numel(ind_sample_tmp) == 0
        continue;
    end
    Nsamplet(id_date,1) = numel(ind_sample_tmp);
    freq_f1_tmp = freq_f1_full(ind_sample_tmp,:);  %Nsample x Nloc
    tc_sub = freq_f1_tmp * pinv(S);  %Nsamplexnica;
    S_sub_tmp = pinv(tc_sub)*freq_f1_tmp;
    S_sub(id_date,:,:) = S_sub_tmp;
    S_sub_tmp_zscore = zeros(size(S_sub_tmp));
    for j = 1:nica
        S_sub_tmp_zscore(j,:) = zscore(S_sub_tmp(j,:));
    end
    S_sub_zscore(id_date,:,:) = S_sub_tmp_zscore;
    save([dual_regression_save_path,'\',date_vec{id_date},'.mat'],...
        'S_sub_tmp_zscore','S_sub_tmp','tc_sub');
end
h1 = figure(123645864); set(h1,'Position',[50,50,2400,600]);
plot(1:Ndate,Nsamplet,'-o','LineWidth',2);
xticks(1:Ndate);xticklabels(date_vec);
title(['Number of Samples per ',num2str(interval_days),'days']);
set(gca,'FontSize',16);
box on; grid on;
saveas(h1,[ica_save_path,'\Nsamplet.jpg']);
end

