function [p_thre,stat_thre] = function_FDR_correction(pval,p_alpha,c,statval)
% controls false positive rate on average to be 0.05;
% sort p-value from the smallest to the largest and find the largest p that
% satisties condition: p(i)<=(i/qmaxM)*alpha/c;
% guassian white noise: c= 1; otherwise give the distribution
    if nargin < 4
        statval = [];
    end
    if nargin < 3
        c = 1;
    end
    [pval_sort,sort_ind] = sort(pval(:),'ascend'); %sort pval in ascending order
    qmaxM = numel(pval);
    seq_compare = (1:qmaxM)'./qmaxM.*p_alpha./c;
    ind_compare = pval_sort<=seq_compare;
    index = find(ind_compare==1);
    if isempty(index)
        p_thre = NaN;   %nothing passed the FDR @ palpha;
        if ~isempty(statval)
            stat_thre = NaN;
        else
            stat_thre = [];
        end
    else
        p_thre = pval_sort(index(end));
        if ~isempty(statval)
            statval_sort = statval(sort_ind);
            stat_thre = statval_sort(index(end));
        else
            stat_thre = [];
        end
    end
end