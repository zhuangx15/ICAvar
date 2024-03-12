clc
clear all
close all

addpath other_Dependence\icasso122\
addpath other_Dependence\FastICA_25\
addpath other_Dependence\common_func\

%% step1
load sample_dataset\covid_test_data.mat freq_f1 freq_f1_name loc_name;  %save data into .mat file with these variables
% freq_f1:  Nsample x Nloc; 
%           alternative allele frequencies at each mutation produced by iVar
% freq_f1_name: Nsample x 3
%               sample-name; sample-collecting-location; sample-collecting-date
% loc_name: Nloc x 1: genome positions
nica = 10; %number of ica component; [your choice];
Niter_icasso = 50; %number of ica-repeat times --> produce stable results
step1_run_ICA(freq_f1,loc_name,nica,Niter_icasso)

%% step2
ica_save_path = [num2str(nica)]; %ica results are saved uncer [current-directory\nica];
freq_f1_date = freq_f1_name(:,3);
iq_threshold = 0.7; %threshold to determine stable ICA components; closer to 1 means more stable
zscore_thre = 2; %zscore-threshold on ICA source file; 2 --> keep 5 positions (following Gaussian distribution);
interval_days = 7; %dual-regression onto weekly samples;
step2_dual_regression(freq_f1,loc_name,freq_f1_date,ica_save_path,iq_threshold,...
    zscore_thre,interval_days);

%% step3
% SARS-CoV-2 Variants of Concerns (VoC) of interest [your choice]; 
variant_of_interest_input = {'B.1.617.2' 'BA.1' 'BA.2' 'BA.2.75' 'BA.2.12.1' 'BA.4' 'BA.5' 'BF.7' 'BQ.1' ...
    'XBB.1' 'XBB.1.5'  'XBB.1.9' 'XBB.1.16' 'XBB.2.3' 'EG.5' 'HV.1' 'FL.1.5.1' 'BA.2.86'};  
% 
barcode_file = 'other_Dependence\prepare_VoC_references_for_annotation\covspectrum_variant_of_interest_hierachical_10182023@0.9_v2.mat';
wuhan_fasta_file =['other_Dependence\wuhan.fasta'];
step3_annotate_dual_regressed_signal(ica_save_path,freq_f1_date,...
    variant_of_interest_input,barcode_file,wuhan_fasta_file,...
    zscore_thre,iq_threshold,interval_days)

%% step4
variant_for_snp_check = {'B.1.617.2','BA.1','BA.2','XBB.1','XBB.1.16','EG.5','HV.1','BA.2.86'}; 
evaluate_SNPs = 'sig-in-gp-ica';
variant_for_snp_path = 'other_Dependence\prepare_VoC_references_for_annotation\covspectrum_dominant_mutations';
step_identify_snp_change_over_time = 1;
step_compute_snp_per_with_significant_change_over_time = 1;
step_limit_potential_novel_pool = 1;

step4_identify_potential_novel_variants(ica_save_path,variant_for_snp_path,variant_for_snp_check,freq_f1_date,...
    wuhan_fasta_file,barcode_file,...
    evaluate_SNPs,interval_days,zscore_thre,iq_threshold,...
    step_identify_snp_change_over_time,step_compute_snp_per_with_significant_change_over_time,step_limit_potential_novel_pool)
