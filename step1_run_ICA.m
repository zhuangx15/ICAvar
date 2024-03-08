function step1_run_ICA(freq_f1,loc_f1,nica,Niter_icasso,root_path)
% Input:
%   freq_f1: NsamplexNlocation; iVar produced alternative allele frequencies at
%   eahc mutation location;
%   loc_f1: Nlocation x 1; mutations;
%   nica: number of ICA component;
%   Niter_icasso: number of times of repeating ICA with different random starting points --> stable results
%   root_path: place to save ica results
% Output:
%   final ica results are saved at: [root_path,'\',num2str(nica)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist("nica",'var')
    nica = 20;
end
if ~exist("Niter_icasso",'var')
    Niter_icasso = 50;
end
if ~exist("root_path",'var')
    root_path = pwd;
end
if ~exist('freq_f1','var') || ~exist("loc_f1",'var')
    error('Please provide sample files');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% remove sequence posistions with max-frequency < Y_max_thre*max across all samples;
ica_save_path = [root_path,'\',num2str(nica)];
if ~exist(ica_save_path,'dir')
    mkdir(ica_save_path);
end
sR=icassoEst('both', freq_f1, Niter_icasso, 'lastEig', nica, 'g', 'pow3', ...
    'approach', 'symm');   %run ICA for Nica_icasso times get stable results;
sR=icassoExp(sR);  %estimate cluster;
[iq,A,W,S,index,h1,h2,h3,h4,h5,h6]=icassoShow(sR,'L', nica,'quality','detailed');
set(h1,"Position",[50,50,1600,1200]);
set(h2,"Position",[50,50,1600,1200]);
set(h3,"Position",[50,50,1600,1200]);
set(h4,"Position",[50,50,1600,1200]);
set(h5,"Position",[50,50,1600,1200]);
set(h6,"Position",[50,50,1600,1200]);
icasso_save_path = [ica_save_path,'\icasso']; 
if ~exist(icasso_save_path,'dir')
    mkdir(icasso_save_path);
end
saveas(h1,[icasso_save_path,'\R_index_determine_number_of_clusters.jpg']);
saveas(h2,[icasso_save_path,'\Estimate_Quality_Iq_for_each_ic.jpg']);
saveas(h3,[icasso_save_path,'\Dendrogram_and_similarity_matrix_among_ic.jpg']);
saveas(h4,[icasso_save_path,'\Similarity_Graph.jpg']);
saveas(h5,[icasso_save_path,'\Source.jpg']);
saveas(h6,[icasso_save_path,'\detailed_estimate_quality.jpg'])
save([ica_save_path,'\ica_results.mat'],'sR','S','A','loc_f1','iq');
end
