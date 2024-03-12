function step4_identify_potential_novel_variants(ica_save_path,variant_for_snp_path,variant_for_snp_check,Yf_sample_date,...
    genome_fasta_file,barcode_file,...
    evaluate_SNPs,interval_days,zscore_thre,ica_iq_threshold,...
    step_identify_snp_change_over_time,step_compute_snp_per_with_significant_change_over_time,step_limit_potential_novel_pool)
%% Input:
%   ica_save_path: directory where ICA results stored;
%   variant_for_snp_path: directory where dominant mutations for each SARS-CoV-2 Variants of concerns (VoC) stored;
%   variant_for_snp_check: VoCs that would check for overlap between their dominant SNPs and ICA-contributing-SNPs
%   Yf_sample_date: sampling date for each wastewater sequencing data
%   barcode_file: reference file containing mutations for each VoC (.mat file required)
%                 .mat could be generated from .txt or .csv file using
%                 function_create_barcode_mat.m
%   genome_fasta_file: wuhan.fasta for SARS-CoV-2 genome; (.fasta file)
%% Input Parameters:
%   evaluate_SNPs: type of SNPs in ica results that are defined as
%   contributing SNPs in ICA (sig-in-gp-ICA: |S|>= zscore) 
%   ica_iq_threshold: threshould to keep stable ICA components (0.7) 
%   zscore_thre: threshold to keep contributing positions in ICA (2)
%   interval_days: intervals to group wastewater samples for dual-regression
%% Input steps:
% step_identify_snp_change_over_time: MAJOR STEP to 
%                                     identify significant SNPs that are changed over time in ICA (set y)
%                                     =1 to GET POTENTIAL NOVEL SNPS                
% step_compute_snp_per_with_significant_change_over_time: PROOF-OF-CONCEPT Step:
%                                                         determine how many in set y is overlapping with dominant mutations in VoCs listed in variant_for_snp_check
% step_limit_potential_novel_pool: PROOF-OF-CONCEPT step:
%                                  determine how many in set y would left after cross-list with known VoCs
%                                  focusing on SNPs with significant
%                                  changes after 'recent_date';
%% Output:
% results saved under:
% [ica_save_path,'\dual_regression_by_',num2str(interval_days),'days\contributing_snps_with_time_effect\',evaluate_SNPs];
%%%%%%%%%%%%%%%%%CAN STILL MODIFY%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist("step_identify_snp_change_over_time",'var')
    % always need to run!!!!
    step_identify_snp_change_over_time = 1;
end
if ~exist('step_compute_snp_per_with_significant_change_over_time','var')
    % this step and the following step are for proof-of-concept
    step_compute_snp_per_with_significant_change_over_time = 1;
end
if ~exist('step_limit_potential_novel_pool','var')  
    step_limit_potential_novel_pool = 1;
end
if step_limit_potential_novel_pool == 1
    % these parameters here are proof-of-concept to show significant
    % time-evolving SNPs identified by ICA after 07/01/2023 are overlapping
    % with 'EG.5,HV.1 or BA.2.86'
    recent_date = datenum('07/01/2023','mm/dd/yyyy');
    variant_keep = {'EG.5','HV.1','BA.2.86'};  %proof of concept
    clustering_method = 'hierachical';
end
addpath other_Dependence\common_func\
if ~exist('variant_for_snp_check','var') || isempty(variant_for_snp_check)
    variant_for_snp_check = {'B.1.617.2','BA.1','BA.2','XBB.1','XBB.1.16','EG.5','HV.1','BA.2.86'}; 
end
if ~exist("zscore_thre",'var')
    zscore_thre = 2;
end
if ~exist('ica_iq_threshold','var')
    ica_iq_threshold = 0.7;
end
if ~exist('interval_days','var')
    interval_days = 7;
end
if ~exist('evaluate_SNPs','var')
    evaluate_SNPs = 'sig-in-gp-ica'; %sig-in-gp-ica; or sig-in-sub-ica or all
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% determine significant SNPs in group ICA map;
load([ica_save_path,'\ica_results.mat'],'S','A','loc_f1','iq');
ind_ica_keep = find(iq>=ica_iq_threshold);  % the estimates are concentrated in m compact and close-to-orthogonal clusters.
% In this case the index to all estimate-clusters is (very close) to one.
% The value drops when the clusters grow wider and mix up
% selected by looking at similarity graph h4
Nsample_in_ica = size(A,1);
S_full = S;
S = S_full(ind_ica_keep,:);
nica = size(S,1);
S_zscore = zeros(size(S));
for j = 1:nica
    S_zscore(j,:) = zscore(S(j,:));
end
S_threshold = abs(S_zscore) >= zscore_thre; %threshold_for_significant_SNPs_on_S;
Nsnp = size(S,2);

%% check overlap between significant SNPs in ICA and SNPs in variant_for_snp_check
N_variant_snp_check = numel(variant_for_snp_check);
snp_check = cell(N_variant_snp_check,1);
for i = 1:numel(variant_for_snp_check)
    T1 = readcell([variant_for_snp_path,'\',variant_for_snp_check{i},'_v2.txt']);
    snp_check{i,1} = T1;
end
snp_detect_f = cell(N_variant_snp_check,1);

%% load dual-regressed data
monthly_type = ['dual_regression_by_',num2str(interval_days),'days'];  %dual_regression_by_month
Yf_sample_date_num = datenum(Yf_sample_date,'mmddyy');
min_date_num = min(Yf_sample_date_num);
max_date_num = max(Yf_sample_date_num);
date_vec_num = min_date_num:interval_days:max_date_num;
date_vec = cellstr(datestr(date_vec_num,'mm-dd-yyyy'));

%% dual regression, similarity to barcode
dual_regression_save_path = [ica_save_path,'\',monthly_type];
Ndate = numel(date_vec);
Nsamplet = nan(Ndate,1);
S_all_date_zscore = nan(Ndate,nica,Nsnp);
for id_date = 1:Ndate
    if ~exist([dual_regression_save_path,'\',date_vec{id_date},'.mat'],'file')
        Nsamplet(id_date,1) = 0;
        continue;
    end
    load([dual_regression_save_path,'\',date_vec{id_date},'.mat'],...
        'S_sub_tmp_zscore','tc_sub');
    S_all_date_zscore(id_date,:,:) = S_sub_tmp_zscore;
    Nsamplet(id_date,1) = size(tc_sub,1);
end
save_path = [dual_regression_save_path,'\contributing_snps_with_time_effect\',evaluate_SNPs];
if ~exist(save_path,'dir')
    mkdir(save_path);
end

ind_non_nan_date = find(Nsamplet>=1);
date_vec_non_nan = date_vec(ind_non_nan_date);
date_vec_num_non_nan = date_vec_num(ind_non_nan_date);
S_non_nan_date_zscore = S_all_date_zscore(ind_non_nan_date,:,:);
Ndate_nonnan = numel(date_vec_num);
ICA_display_name = strcat('ICA-', strsplit(num2str(1:nica),' '));
t_Vec = date_vec_num_non_nan - min(date_vec_num_non_nan);
t_Vec = t_Vec(:);

if step_identify_snp_change_over_time == 1
    flag_add_constant = 1;
    sig_snp_with_sig_time_change_name_f = [];
    sig_snp_name_f = [];
    snp_name_sig_change_over_time = cell(nica,1);
    pthre_method_t = cell(nica,1);
    S_all_sub_tmp_sig_t = cell(nica,1);
    S_snp_recode_del_t = cell(nica,3);
    for i = 1:nica
        switch evaluate_SNPs
            case 'sig-in-gp-ica'
                ind_sig_SNP = find(S_threshold(i,:)==1);
            case 'sig-in-sub-ica'
                all_sub_this_ica = squeeze(S_non_nan_date_zscore(:,i,:));
                ind_sig_SNP = find(max(abs(all_sub_this_ica))>=zscore_thre);
            case 'all'
                ind_sig_SNP = 1:Nsnp;
        end
        snp_name = loc_f1(ind_sig_SNP);
        sig_snp_name_f = [sig_snp_name_f;snp_name(:)];

        S_all_sub_tmp = abs(squeeze(S_non_nan_date_zscore(:,i,ind_sig_SNP)));
        [R_square_model,Bfull,regressor_T,regressor_P,regressor_D,...
            con_con_map,con_t_map,con_p_map,con_d_map,rho_map,var_con_con_map,...
            RES,dfe,Xfull,SSE,SSR,SST,F_model,p_model,ind_empty_col_in_X] = XZ_SPM_univariate_analysis...
            (t_Vec,S_all_sub_tmp,[],flag_add_constant);
        [fdr_thre] = function_FDR_correction(regressor_P(2,:),0.05);
        fwe_thre = 0.05/numel(snp_name);
        if isnan(fdr_thre)
            pthre = 0.05;
            pthre_method = 'raw';
        elseif fdr_thre>fwe_thre
            pthre = fwe_thre;
            pthre_method = 'fwe';
        else
            pthre = fdr_thre;
            pthre_method = 'fdr';
        end

        ind = find(regressor_P(2,:)<=pthre);
        S_all_sub_tmp_sig = S_all_sub_tmp(:,ind);
        snp_name_sig_p = snp_name(ind);
        sig_p = regressor_P(2,ind);
        snp_name_sig_change_over_time{i,1} = snp_name_sig_p;
        snp_name_sig_change_over_time{i,2} = sig_p;
        S_all_sub_tmp_sig_t{i,1} = S_all_sub_tmp_sig;
        sig_snp_with_sig_time_change_name_f = [sig_snp_with_sig_time_change_name_f;snp_name_sig_p(:)];

        pthre_method_t{i,1} = pthre_method;
        if strcmp(pthre_method,'raw') && numel(ind)>=10
            disp(['ica-',num2str(i),'not-survive-correction']);
            continue;
        end
        h988 = figure(1560);
        Nsig = size(S_all_sub_tmp_sig,2);
        if Nsig >= 20
            height = 1200;
            fontsize_tmp = 12;
        elseif Nsig >= 10
            height = 800;
            fontsize_tmp = 16;
        else
            height = 600;
            fontsize_tmp = 16;
        end
        set(h988,'Position',[50,50,3200,height]);
        h998 = signalplot(S_all_sub_tmp_sig',h988);
        xticks(1:Ndate_nonnan); xticklabels(date_vec_non_nan); box on; grid on;
        yticks(1:numel(ind));yticklabels(snp_name_sig_p);
        title(sprintf([pthre_method,'-corrected; raw-p<%.02e'],pthre));
        set(gca,'FontSize',fontsize_tmp);
        saveas(h998,[save_path,'\ica',num2str(i),'_threshold_significant_change_over_time.jpg']);

        %%recode deletion
        [sig_p_recode_del,snp_name_sig_p_recode_del] = ...
            recode_deletion(genome_fasta_file,snp_name_sig_p,sig_p);
        [S_all_sub_tmp_sig_recode_del,snp_name_sig_p_recode_del1] = ...
            recode_deletion(genome_fasta_file,snp_name_sig_p,S_all_sub_tmp_sig);
        sum(~strcmp(snp_name_sig_p_recode_del,snp_name_sig_p_recode_del1));
        S_snp_recode_del_t{i,1} = snp_name_sig_p_recode_del;  %snp name after recode
        S_snp_recode_del_t{i,2} = sig_p_recode_del; %pvalue at recode snp name
        S_snp_recode_del_t{i,3} = S_all_sub_tmp_sig_recode_del; %S for plot
    end


    for ii = 1:N_variant_snp_check
        tmp_snp_check = snp_check{ii};
        snp_detect = nan(nica,numel(tmp_snp_check)); %put pvalue here
        for i = 1:nica
            snp_name_sig_p_recode_del = S_snp_recode_del_t{i,1};
            sig_p_recode_del = S_snp_recode_del_t{i,2};
            for jj = 1:numel(tmp_snp_check)
                ind_detect = find(contains(snp_name_sig_p_recode_del,tmp_snp_check{jj}));
                if numel(ind_detect)==1
                    snp_detect(i,jj) = sig_p_recode_del(ind_detect);
                end
            end
        end
        snp_detect_f{ii,1} = snp_detect;
        xlsout = [[{''};ICA_display_name(:)] [tmp_snp_check(:)';num2cell(snp_detect)]];
        % SNPs in checked VoC detected by ICA
        xlswrite([save_path,'\detection_matrix.xlsx'],xlsout,variant_for_snp_check{ii});

        %% plot
        [min_tmp_detect,ind_ica_min_p] = min(snp_detect);
        ind_detect_snp = find(~isnan(min_tmp_detect));
        ind_ica_min_p = ind_ica_min_p(ind_detect_snp);
        S_plot = zeros(numel(t_Vec),numel(ind_detect_snp));
        for jj = 1:numel(ind_detect_snp)
            current_snp = tmp_snp_check{ind_detect_snp(jj)};
            ind = find(contains(S_snp_recode_del_t{ind_ica_min_p(jj),1},current_snp));
            S_plot(:,jj) = S_snp_recode_del_t{ind_ica_min_p(jj),3}(:,ind);
        end
        Nsig = numel(ind_detect_snp);
        h9881 = figure(156000);
        if Nsig >= 20
            height = 1200;
            fontsize_tmp = 16;
        elseif Nsig >= 10
            height = 800;
            fontsize_tmp = 16;
        else
            height = 600;
            fontsize_tmp = 16;
        end
        set(h9881,'Position',[50,50,3200,height]);
        h9981 = signalplot(S_plot',h9881);
        xticks(1:Ndate_nonnan); xticklabels(date_vec_non_nan); box on; grid on;
        yticks(1:Nsig);yticklabels(tmp_snp_check(ind_detect_snp));
        title(variant_for_snp_check{ii});
        set(gca,'FontSize',fontsize_tmp);
        set(gca,'FontWeight','bold');
        saveas(h9981,[save_path,'\min_sig_p_from_all_ica_significant_change_over_time_',variant_for_snp_check{ii},'.jpg']);
    end
    sig_snp_name_f = unique(sig_snp_name_f);
    sig_snp_with_sig_time_change_name_f = unique(sig_snp_with_sig_time_change_name_f);
    save([save_path,'\sig_SNPs_over_time.mat'],'snp_name_sig_change_over_time','S_all_sub_tmp_sig_t','pthre_method_t',...
        'variant_for_snp_check','snp_check','S_snp_recode_del_t',...
        'sig_snp_with_sig_time_change_name_f','sig_snp_name_f');
end

if step_compute_snp_per_with_significant_change_over_time == 1
    xlsout_detect_table = zeros(numel(variant_for_snp_check),3);
    for ii = 1:numel(variant_for_snp_check)
        [~,~,RAW] = xlsread([save_path,'\detection_matrix.xlsx'],variant_for_snp_check{ii});
        snp = RAW(1,2:end);
        Nsnp = numel(snp);
        RAW_detect = cell2mat(RAW(2:end,2:end));
        ind = find(sum(isnan(RAW_detect))<size(RAW_detect,1));
        Ndetect = numel(ind);
        xlsout_detect_table(ii,:) = [Nsnp,Ndetect,Ndetect/Nsnp];
    end
    xlsout_detect_table_f = [{'Variants','Total-snp','detected-snp','detected-percentage'};[variant_for_snp_check(:) num2cell(xlsout_detect_table)]];
    xlswrite([save_path,'\detection_matrix.xlsx'],xlsout_detect_table_f,'final');
end


if step_limit_potential_novel_pool == 1  %clustering
    %% potential novel: first eliminate known ones;
%     clustering_method = 'hierachical';
    load([save_path,'\sig_SNPs_over_time.mat'],'sig_snp_with_sig_time_change_name_f','S_snp_recode_del_t');
    potential_novel_save_path = [save_path,'\potential_novel']; 
    if ~exist(potential_novel_save_path,'dir')
        mkdir(potential_novel_save_path);
    end
    [~,sig_snp_with_sig_time_change_name_f_recode_del,feature_remove] = recode_deletion...
        (genome_fasta_file,sig_snp_with_sig_time_change_name_f);
    % eliminate known ones;
    load(barcode_file,'loc_barcode');
    % but keep the ones needed;
%     variant_keep = {'EG.5','HV.1','BA.2.86'};  %proof of concept
    [~,ind11,ind22] = intersect(variant_for_snp_check,variant_keep,'stable');
    snp_variant = snp_check(ind11);
    snp_variant_f = [];
    for j = 1:numel(variant_keep)
        snp_variant_f = [snp_variant_f;snp_variant{j}];
    end
    snp_variant_f = unique(snp_variant_f);
    [loc_barcode_known] = setdiff(loc_barcode,snp_variant_f);
    [potential_novel,ind1] = setdiff(sig_snp_with_sig_time_change_name_f_recode_del,loc_barcode_known);

%     recent_date = datenum('07/01/2023','mm/dd/yyyy');
    ind_recent = find(date_vec_num_non_nan>=recent_date);
    t_Vec_recent = t_Vec(ind_recent);

    S_potential_novel_all_ICA = nan(numel(t_Vec),numel(potential_novel));  %Nt x Nsnp --> Nfeature x Nsample
    S_potential_novel_all_ICA_zscore = nan(numel(t_Vec),numel(potential_novel));
    for j = 1:numel(potential_novel)
        S_tmp_all_ica = zeros(numel(t_Vec),nica);
        S_tmp_all_ica_zscore = zeros(numel(t_Vec),nica);

        snp_p_tmp_all_ica = ones(nica,1);
        for i = 1:nica
            snp_name_tmp = S_snp_recode_del_t{i,1};
            snp_p_tmp = S_snp_recode_del_t{i,2};
            S_tmp = S_snp_recode_del_t{i,3};
            ind_this_snp = find(strcmp(snp_name_tmp,potential_novel{j}));
            if numel(ind_this_snp)==1
                S_tmp_all_ica(:,i) = S_tmp(:,ind_this_snp);
                S_tmp_all_ica_zscore(:,i) = zscore(S_tmp(:,ind_this_snp));
                snp_p_tmp_all_ica(i,1) = snp_p_tmp(ind_this_snp);
            elseif numel(ind_this_snp)>1
                disp(['MoreThanOne2@',potential_novel{j}]);
                [~,ind] = min(snp_p_tmp(ind_this_snp));
                S_tmp_all_ica(:,i) = S_tmp(:,ind_this_snp(ind));
                snp_p_tmp_all_ica(i,1) = snp_p_tmp(ind_this_snp(ind));
                S_tmp_all_ica_zscore(:,i)  = zscore(S_tmp(:,ind_this_snp(ind)));
            end
        end
        [~,ind_min] = min(snp_p_tmp_all_ica);
        S_potential_novel_all_ICA(:,j) = S_tmp_all_ica(:,ind_min);  %18 tp x 113 snp;
        S_potential_novel_all_ICA_zscore(:,j) = S_tmp_all_ica_zscore(:,ind_min);
    end
    S_potential_novel_all_ICA_zscore = S_potential_novel_all_ICA_zscore(ind_recent,:);
    if ~exist([potential_novel_save_path,'\',clustering_method],'dir')
        mkdir([potential_novel_save_path,'\',clustering_method]);
    end
    switch clustering_method
        case 'hierachical'
            %% hierachical clustering
            Z = linkage(S_potential_novel_all_ICA_zscore','ward');  %after normalization should use Euclidean distance
            h3 = figure (10096); set(h3,'Position',[50,50,1600,2400]);
            dendrogram_oritation = 'left';
            [~,~,OUTPERM] = dendrogram(Z,0,'Orientation',dendrogram_oritation,'Labels',potential_novel); xtickangle(90);
            box on; grid on;
            saveas(h3,[potential_novel_save_path,'\',clustering_method,'\dendrogram.jpg']);

            eva = evalclusters(S_potential_novel_all_ICA_zscore','linkage','silhouette','KList',[2:10]); %euclidean distance; already zscored;;
            %             Ncluster = eva.OptimalK;
            %             Ncluster = 7;
            T = cluster(Z,'cutoff',18,'criterion','distance');
            Ncluster = numel(unique(T));
            %             T = cluster(Z,"maxclust",Ncluster);
        case 'kmeans'
            eva = evalclusters(S_potential_novel_all_ICA_zscore','kmeans','gap','KList',[2:20]);
    end
    save_path_f = [potential_novel_save_path,'\',clustering_method,'\',num2str(Ncluster),'Clusters'];
    if ~exist(save_path_f,'dir')
        mkdir(save_path_f);
    end

    T_snp = cell(Ncluster,1);

    overlap_with_variant_keep = zeros(numel(variant_keep),Ncluster);
    overlap_with_variant_keep_detail = cell(numel(variant_keep),Ncluster);
    for i = 1:Ncluster
        % plot
        S_tmp = S_potential_novel_all_ICA(:,T==i);
        T_snp{i,1} = potential_novel(T==i);
        Ntmp = numel(T_snp{i,1});
        if Ntmp>=10
            height = 1200;
        else
            height = 600;
        end
        h98811 = figure(500123);
        set(h98811,'Position',[50,50,3200,height]);
        h9981 = signalplot(S_tmp',h98811);

        xticks(1:Ndate_nonnan); xticklabels(date_vec_non_nan); box on; grid on;
        yticks(1:numel(T_snp{i,1}));yticklabels(T_snp{i,1});
        title(['Cluster-',num2str(i)]);
        set(gca,'FontSize',22);
        set(gca,'FontWeight','bold');
        saveas(h98811,[save_path_f,'\Cluster',num2str(i),'.jpg']);

        for j = 1:numel(variant_keep)
            [overlap_snp,ind11,ind22] = intersect(snp_variant{j},T_snp{i,1});
            overlap_with_variant_keep(j,i) = numel(overlap_snp);
            overlap_with_variant_keep_detail{j,i} = overlap_snp;
        end
    end
    xlsout = [[{''},strcat('cluster-', strsplit(num2str(1:Ncluster),' '))];[variant_keep(:) num2cell(overlap_with_variant_keep)]];
    save([save_path_f,'\overlap_with_known_variants.mat'],'xlsout','T_snp','T','overlap_with_variant_keep_detail');
end
disp(['Results saved under: ',save_path]);
end