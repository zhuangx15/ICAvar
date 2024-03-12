function step3_annotate_dual_regressed_signal(ica_save_path,Yf_sample_date,...
    variant_of_interest_input,barcode_file,wuhan_fasta_file,...
    zscore_thre,ica_iq_threshold,interval_days)
%% Input:
%   ica_save_path: directory where ICA results stored;
%   Yf_sample_date: sampling date for each wastewater sequencing data
%   variant_of_interest_input: SARS-CoV-2 VoCs (just name);
%   barcode_file: reference file containing mutations for each VoC (.mat file required)
%                 .mat could be generated from .txt or .csv file using
%                 function_create_barcode_mat.m
%   wuhan_fasta_file: wuhan.fasta for SARS-CoV-2 genome; (.fasta file)
%% Input Parameters:
%   ica_iq_threshold: threshould to keep stable ICA components (0.7) 
%   zscore_thre: threshold to keep contributing positions in ICA (2)
%   interval_days: intervals to group wastewater samples for dual-regression
%% Output:
% results saved under:
% [ica_save_path,'\dual_regression_by_',num2str(interval_days),'days\detect_matrix_ICA'];
%% Additional parameters you could modify here:
%%=================threshold to determine whether a VoC is detected or not
%%follows the description in Zhuang et al., 2024, but you could still modify them===
    level1_sen_thre = 0.5;  %barcode_level == 1 & Nsnp_in_each_variant>=100
    level1dot5_F1_thre = 0.3; %barcode_level >1 & Nsnp_in_each_variant > 20
    level2_sen_thre = 0.6; %barcode_level >1 & Nsnp_in_each_variant: [10,20];
    level3_sen_thre = 0.5; %barcode_level >1 & Nsnp_in_each_variant: (0,10);
% ===================================================================
    if ~exist("zscore_thre",'var')
        zscore_thre = 2;
    end
    if ~exist('ica_iq_threshold','var')
        ica_iq_threshold = 0.7;
    end
    if ~exist('interval_days','var')
        interval_days = 7; 
    end
    addpath other_Dependence\common_func; %for recode.m; compute_binary_AUC.m     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%d
    close all
    load(barcode_file,"loc_barcode","barcode_f","barcode_level","barcode_parent","Variant_of_interest");
    barcode_variant=Variant_of_interest; clear "Variant_of_interest";
    disp(['current-ICAvar-save-path: ',ica_save_path]);
    monthly_type = ['dual_regression_by_',num2str(interval_days),'days'];  %dual_regression_by_month
    dual_regression_save_path = [ica_save_path,'\',monthly_type];
    final_detection_save_path = [dual_regression_save_path,'\detect_matrix_ICA'];  %detect_matrix_ICA_v2_add_level4

    if ~exist(final_detection_save_path,'dir')
        mkdir(final_detection_save_path);
    end
    %% dual regression already done; load results; dual regress every day
    % for variant order;
    [Variant_of_interest_ICA,ind1,ind2] = intersect(variant_of_interest_input,barcode_variant,'stable');
    barcode_variant = barcode_variant(ind2);  %recorder barcode to follow VoC-input;
    barcode_f = barcode_f(ind2,:);
    barcode_level = barcode_level(ind2,:);
    barcode_parent = barcode_parent(ind2,:);
    Nsnp = numel(loc_barcode);
    
    Nvariant_ICA = numel(Variant_of_interest_ICA);
    load([ica_save_path,'\ica_results.mat'],'S','A','loc_f1','iq');
    ind_ica_keep = find(iq>=ica_iq_threshold);  % the estimates are concentrated in m compact and close-to-orthogonal clusters.
                                                % In this case the index to all estimate-clusters is (very close) to one.
                                                % The value drops when the clusters grow wider and mix up
                                                % selected by looking at similarity graph h4
                                                % recode deletion
    [S,loc_f1_recode_del] = recode_deletion(wuhan_fasta_file,loc_f1,S);
    Nsample_in_ica = size(A,1);
    S_full = S;
    S = S_full(ind_ica_keep,:);
    nica = size(S,1);
    S_zscore = zeros(size(S));
    for j = 1:nica
        S_zscore(j,:) = zscore(S(j,:));
    end
    S_threshold = abs(S_zscore) >= zscore_thre; %threshold_for_significant_SNPs_on_S;
    
    %% intersect of barcode and sample location
    [~,ind_barcode,ind_f1] = intersect(loc_barcode,loc_f1_recode_del,'stable');
    barcode_f = barcode_f(:,ind_barcode);
    loc_barcode = loc_barcode(ind_barcode);
    disp([num2str(numel(ind_barcode)),'/',num2str(Nsnp),'-position-in-barcode-present-in-samples']);
    
    S_f2_threshold = double(S_threshold(:,ind_f1));
    S_f2 = S(:,ind_f1);
    S_zscore_f2 = S_zscore(:,ind_f1);
    Yf_sample_date_num = datenum(Yf_sample_date,'mmddyy');
    min_date_num = min(Yf_sample_date_num);
    max_date_num = max(Yf_sample_date_num);
    min_date = datestr(min_date_num,'mmddyy');
    max_date = datestr(max_date_num,'mmddyy');
    date_vec_num = min_date_num:interval_days:max_date_num;
    date_vec = cellstr(datestr(date_vec_num,'mm-dd-yyyy'));
    %% Annotation: similarity to barcode
    Ndate = numel(date_vec);
    Nsamplet = nan(Ndate,1);
    detect_matrix_ICA_corr = nan(Nvariant_ICA,Ndate);
    detect_matrix_ICA_corr_p = nan(Nvariant_ICA,Ndate);
    detect_matrix_ICA_JI = nan(Nvariant_ICA,Ndate);
    detect_matrix_ICA_sensitivity = nan(Nvariant_ICA,Ndate);
    detect_matrix_ICA_specificity = nan(Nvariant_ICA,Ndate);
    detect_matrix_ICA_AUC = nan(Nvariant_ICA,Ndate);
    detect_matrix_ICA_F1score = nan(Nvariant_ICA,Ndate);
    for id_date = 1:Ndate
        if ~exist([dual_regression_save_path,'\',date_vec{id_date},'.mat'],'file')
            Nsamplet(id_date,1) = 0;
            continue;
        end
        load([dual_regression_save_path,'\',date_vec{id_date},'.mat'],...
            'S_sub_tmp_zscore','tc_sub');
        Nsamplet(id_date,1) = size(tc_sub,1);
        %%  hierachical similarity to barcode
        % recode deletion in ICA
        [S_sub_tmp_zscore,loc_f1_recoded1] = recode_deletion(wuhan_fasta_file,loc_f1,S_sub_tmp_zscore);
        if ~sum(strcmp(loc_f1_recode_del,loc_f1_recoded1)) == numel(loc_f1_recoded1)
            warning([date_vec{id_date},' position not match group map']);
        end
        S_sub_tmp_f2 = S_sub_tmp_zscore(:,ind_f1);
        S_sub_tmp_f2_binary = double(abs(S_sub_tmp_f2)>=2);
        % compare against barcode_f
        [AUC,sensitivity,specificity,F1score] = compute_AUC_binary(S_sub_tmp_f2_binary,barcode_f);  %nica x Nvariant;
        [this_date_ICA_detect_sensitivity,ind_max_sen] = max(sensitivity); %1xNvariant
        this_date_ICA_detect_AUC = nan(Nvariant_ICA,1);
        this_date_ICA_detect_specificity = nan(Nvariant_ICA,1);
        this_date_ICA_detect_F1score = nan(Nvariant_ICA,1);
        for ii = 1:Nvariant_ICA
            this_date_ICA_detect_AUC(ii,1) = AUC(ind_max_sen(ii),ii); %1xNvariant
            this_date_ICA_detect_specificity(ii,1) = specificity(ind_max_sen(ii),ii); %1xNvariant
            this_date_ICA_detect_F1score(ii,1) = F1score(ind_max_sen(ii),ii); %1xNvariant
        end
        [JI_sub,corr_sub,p_corr_sub,ind_loc_overlap] = compute_JI_and_corr(S_sub_tmp_f2_binary,barcode_f);  %nica x Nvariant;
        this_date_ICA_detect_JI = max(JI_sub); %1xNvariant
        [this_date_ICA_detect_corr,ind_max_corr] = max(corr_sub); %1xNvariant
        this_date_ICA_detect_corr_p = nan(Nvariant_ICA,1);
        for ii = 1:Nvariant_ICA
            this_date_ICA_detect_corr_p(ii,1) = p_corr_sub(ind_max_corr(ii),ii); %1xNvariant
        end
        detect_matrix_ICA_AUC(:,id_date) = this_date_ICA_detect_AUC;
        detect_matrix_ICA_corr_p(:,id_date) = this_date_ICA_detect_corr_p;
        detect_matrix_ICA_corr(:,id_date) = this_date_ICA_detect_corr;
        detect_matrix_ICA_JI(:,id_date) = this_date_ICA_detect_JI;
        detect_matrix_ICA_sensitivity(:,id_date) = this_date_ICA_detect_sensitivity;
        detect_matrix_ICA_specificity(:,id_date) = this_date_ICA_detect_specificity;
        detect_matrix_ICA_F1score(:,id_date) = this_date_ICA_detect_F1score;
    end
    save([final_detection_save_path,'\detect_matrix_ICA_hierachical.mat'],'date_vec','Nsamplet',...
        'detect_matrix_ICA_JI','detect_matrix_ICA_corr','detect_matrix_ICA_corr_p',...
        'detect_matrix_ICA_sensitivity','detect_matrix_ICA_F1score','detect_matrix_ICA_specificity','detect_matrix_ICA_AUC',...
        'Variant_of_interest_ICA','barcode_f','loc_barcode','barcode_parent','barcode_level');


    %% determine detection matrix: thresholding annotation matrices
    detect_matrix_name = {'Spearman-corr','JI','sensitivity','specificity','F1score','AUC'};
    Neval = numel(detect_matrix_name);
    h1 = figure(15456);
    scatter(detect_matrix_ICA_JI(:),detect_matrix_ICA_F1score(:));
    box on; grid on;
    xlabel('JI'); ylabel('F1score');
    set(gca,'FontSize',16);
    saveas(h1,[final_detection_save_path,'\JI_vs_F1score.jpg']);

    detect_matrix = zeros(Neval,Nvariant_ICA,Ndate);
    detect_matrix(1,:,:) = detect_matrix_ICA_corr;
    detect_matrix(2,:,:) = detect_matrix_ICA_JI;
    detect_matrix(3,:,:) = detect_matrix_ICA_sensitivity;
    detect_matrix(4,:,:) = detect_matrix_ICA_specificity;
    detect_matrix(5,:,:) = detect_matrix_ICA_F1score;
    detect_matrix(6,:,:) = detect_matrix_ICA_AUC;

    ind_date_display = find(Nsamplet>=1);
    for i = 1:Neval
        [h32] = plot_detecting_matrix(squeeze(detect_matrix(i,:,:)),[],date_vec,ind_date_display,...
            Variant_of_interest_ICA,[],detect_matrix_name{i});
        saveas(h32,[final_detection_save_path,'\detected_matrix_ICA_',detect_matrix_name{i},'.jpg']);
    end
    detect_matrix_q = reshape(detect_matrix,[Neval,Nvariant_ICA*Ndate]);
    ind_nonNan = find(sum(isnan(detect_matrix_q))==0);
    detect_matrix_q = detect_matrix_q(:,ind_nonNan);

    Y = pdist(detect_matrix_q,"correlation");
    Z = linkage(Y);
    h10086 = figure(45875213); dendrogram(Z,0,'labels',detect_matrix_name); box on; grid on;
    saveas(h10086,[final_detection_save_path,'\similarity_between_evalution_matrix.jpg']);
   
    fwe_p = 0.05/numel(detect_matrix_ICA_corr_p);
    detect_matrix_ICA_corr_p_fdr = mafdr(detect_matrix_ICA_corr_p(:),'BHFDR',true);
    detect_matrix_ICA_corr_p_fdr = reshape(detect_matrix_ICA_corr_p_fdr,[Nvariant_ICA,Ndate]);


    min_corr_pass_fwe_p = min(abs(detect_matrix_ICA_corr(detect_matrix_ICA_corr_p<=fwe_p)));
    detect_matrix_f = zeros(size(detect_matrix_ICA_sensitivity));
    Nsnp_in_each_variant = sum(barcode_f,2);

    % level1: barcode_level == 1
    ind_level11 = find(barcode_level==1 & Nsnp_in_each_variant<100);
    detect_matrix_level11_corr_p = detect_matrix_ICA_corr_p(ind_level11,:);
    detect_matrix_level11 = (detect_matrix_level11_corr_p<=fwe_p);
    detect_matrix_f(ind_level11,:) = detect_matrix_level11;

    ind_level12 = find(barcode_level==1 & Nsnp_in_each_variant>=100);
    detect_matrix_level12_corr_p = detect_matrix_ICA_corr_p(ind_level12,:);
    detect_matrix_level12_sensitivity = detect_matrix_ICA_sensitivity(ind_level12,:);
    detect_matrix_level12 = (detect_matrix_level12_corr_p<=fwe_p & detect_matrix_level12_sensitivity>level1_sen_thre);  %BA.1
    detect_matrix_f(ind_level12,:) = detect_matrix_level12;
    level1_criteria = ['(level1Variant&Nsnp<100-with-fweCorr)|(level1Variant&Nsnp>=200)-with-fweCorr&senitivity>',num2str(level1_sen_thre)];

    % level1dot5: variants with determined mutations > 20
    ind_level1dot5 = find(Nsnp_in_each_variant>20 & barcode_level>1);  
    detect_matrix_level1dot5_corr_p = detect_matrix_ICA_corr_p(ind_level1dot5,:);
    detect_matrix_level1dot5_corr_p_fdr = detect_matrix_ICA_corr_p_fdr(ind_level1dot5,:);
    detect_matrix_level1dot5_F1score = detect_matrix_ICA_F1score(ind_level1dot5,:);
    detect_matrix_level1dot5 = detect_matrix_level1dot5_corr_p_fdr<=0.05 & detect_matrix_level1dot5_F1score>=level1dot5_F1_thre;
    detect_matrix_f(ind_level1dot5,:) = detect_matrix_level1dot5;
    level1dot5_criteria = ['level234Variants&moreThan20mutation-with-fdrCorr&F1score>=',num2str(level1dot5_F1_thre)];

    % level2: [10,20]; sensitivity >= 0.6;
    ind_level2 = find(Nsnp_in_each_variant<=20 & Nsnp_in_each_variant>=10);
    detect_matrix_level2_sensitivity = detect_matrix_ICA_sensitivity(ind_level2,:);
    detect_matrix_level2  = detect_matrix_level2_sensitivity>=level2_sen_thre;
    detect_matrix_f(ind_level2,:) = detect_matrix_level2;
    level2_criteria = ['level234Variants&[10,20]mutation-with-sensitivity>=',num2str(level2_sen_thre)];

    % level3: (0,10); sensitivity >= 0.5
    ind_level3 = find(Nsnp_in_each_variant<10);
    detect_matrix_level3_sensitivity = detect_matrix_ICA_sensitivity(ind_level3,:);
    detect_matrix_level3  = detect_matrix_level3_sensitivity>=level3_sen_thre;
    detect_matrix_f(ind_level3,:) = detect_matrix_level3;
    level3_criteria = ['level234Variants&(0,10)mutation-with-sensitivity>=',num2str(level3_sen_thre)];
    criteria = {level1_criteria;level1dot5_criteria;level2_criteria;level3_criteria};

    %% detect matrix
    [h3,xlsout_first_last_date] = plot_detecting_matrix([],detect_matrix_f,date_vec,ind_date_display,...
        Variant_of_interest_ICA,[],'ICA-final-threshold');
    [Nr,Nc] = size(xlsout_first_last_date);
    xlsout_first_last_date{Nr+1,1} = 'criteria';
    xlsout_first_last_date(Nr+2:Nr+numel(criteria)+1,1) = criteria;
    saveas(h3,[final_detection_save_path,'\detected_by_ICA_hierachical_threshold.jpg']);
    xlswrite([final_detection_save_path,'\first_last_date_new_criteria_11132023.xlsx'],xlsout_first_last_date);
    xlsout_detect_matrix_binary = [[{''};Variant_of_interest_ICA(:)] [date_vec(:)';num2cell(detect_matrix_f)]];
    xlswrite([final_detection_save_path,'\detect_matrix_f.xlsx'],xlsout_detect_matrix_binary,'binary');
    %% remove dates with only 1 sample
    ind_more_than_1_sample = find(Nsamplet>1);
    [h3,xlsout_first_last_date] = plot_detecting_matrix([],detect_matrix_f,date_vec,ind_more_than_1_sample,...
        Variant_of_interest_ICA,[],'ICA-final-threshold');
    saveas(h3,[final_detection_save_path,'\detected_by_ICA_hierachical_threshold_remove_Nsamplet1.jpg']);
    xlswrite([final_detection_save_path,'\first_last_date_new_criteria_11132023.xlsx'],xlsout_first_last_date,'remove_Nsamplet=1');
    disp(['Annotated results saved under: ', final_detection_save_path]);
end
