# ICAvar
Early Detection of Novel SARS-CoV-2 Variants from Urban and Rural Wastewater through Genome Sequencing and Machine Learning 

## Software and Packages required: 
1. Matlab: https://www.mathworks.com/
2. ICAsso*: https://research.ics.aalto.fi/ica/icasso/
3. FastICA*: https://research.ics.aalto.fi/ica/fastica/ <br>
 *Under other_Dependence folder we have provided icasso122 and FastICA_25 packages for the purpose of testing the code. You MUST obtain proper usage permissions from the orignial authors. 

## Steps:
1. Running ICA
2. Dual regression
3. Annotate dual-regressed signal to known COVID strains
4. Identify potential novel mutations

## Testing dataset:
script_example_code_test_all_steps.m provides an example of running all four steps with sample_dataset\covid_test_data.mat <br>
Results reported in Zhuang et al., 2024 were produced using sample_dataset\Yf_all_ivar_variant_sep21_nov23_50xcoverage_gt_80_vRefine_vNoDup_github.mat
