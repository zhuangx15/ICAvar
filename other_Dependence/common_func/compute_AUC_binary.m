function [AUC,sensitivity,specificity,F1score] = compute_AUC_binary(S_f2_threshold,X1_inclusive)
Nica = size(S_f2_threshold,1);
Nvariant_in_barcode = size(X1_inclusive,1);
AUC = nan(Nica, Nvariant_in_barcode);
sensitivity = nan(Nica,Nvariant_in_barcode);
specificity = nan(Nica,Nvariant_in_barcode);
F1score = nan(Nica,Nvariant_in_barcode);
for i = 1:Nica
    for j = 1:Nvariant_in_barcode
        s = scoreCal(S_f2_threshold(i,:)',X1_inclusive(j,:)');
        AUC(i,j) = s.AUC;
        sensitivity(i,j) = s.sensitivityORrecall;
        specificity(i,j) = s.specificity;
        F1score(i,j) = s.F1score;
    end
end
end


function [s1,jaccard] = scoreCal(Y_vali_f1,Y,flag_plot)
if nargin < 3
    flag_plot = 0;
end
Y_vali_f1(Y_vali_f1~=1) = 0;
Y(Y~=1)=0;
N0 = size(Y,1);

TP = sum(Y_vali_f1==1 & Y==1);
FP = sum(Y_vali_f1==1 & Y==0);
TN = sum(Y_vali_f1==0 & Y==0);
FN = sum(Y_vali_f1==0 & Y==1);
FPR = FP/(FP+TN+eps);    %false positive rate
FNR = FN/(TP+FN+eps);
acc = (TP+TN) / N0;
precision = TP./(TP+FP);
sensitivity = TP/(TP+FN);  %also called recall
%F1_score = 2*(sensitivity*precision)/(sensitivity+precision);
F1_score = (2*TP)/(2*TP+FP+FN);
jaccard = TP / (TP + FP + FN);
FNplusFP = FPR + FNR;
specificity = TN/(TN+FP); 
if flag_plot == 1
    [~,~,~,AUC1] = perfcurve(Y,double(Y_vali_f1),1);
else
    AUC1 = (1-specificity)*sensitivity/2 + (sensitivity+1)*specificity/2;
end
s1 = struct('acc',acc,'precision',precision,'sensitivityORrecall',sensitivity,'specificity',specificity,...
    'F1score',F1_score,'AUC',AUC1,'FNandFP',FNplusFP);
end