clc
clear all
close all
family_tree_file = 'family_tree_barcode_of_interest_delta_omicron_github.xlsx';
cov_spectrum_path = pwd;
covspectrum_thre = 0.9;
prepare_barcode_covspectrum_v2(family_tree_file,cov_spectrum_path,covspectrum_thre)