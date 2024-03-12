function [Jaccard_Index_level1,corr_Index_level1,p_corr,ind] = compute_JI_and_corr(S_f2_threshold,X1_inclusive)
A1 = [S_f2_threshold;X1_inclusive];
ind = find(sum(A1==0)<size(A1,1));
Jaccard_Index_level1 = 1-pdist2(S_f2_threshold(:,ind),X1_inclusive(:,ind),'jaccard');
[corr_Index_level1,p_corr] = corr(S_f2_threshold(:,ind)',X1_inclusive(:,ind)','type','Spearman');
end