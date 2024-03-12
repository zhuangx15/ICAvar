function pval = tval_2_pval(tval,df,tail)
switch tail
    case 'right'
        pval = 1 - tcdf(tval,df);
    case 'both'
        pval = 2 * tcdf(-abs(tval), df);
    case 'left'
        pval = tcdf(tval, df);
end
% pval(pval==0) = eps;
end