function [X1_inclusive,variant_level1_inclusive,X1_only,variant_level1_only] = determine_barcode_for_each_hierachical_level(RAW,barcode,variant_name)
% RAW is variants (with its sublineages) for this level;
% this-level variant are the first element before ',' in RAW
% sublineages are elements after ','; %use strsplit for this purpose
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nr = size(RAW,1);
Nsnp = size(barcode,2);
level1_variant = cell(Nr,1);
X1tmp = zeros(Nr,Nsnp);  %X for all variants under covspectrum; sublineages also has one row; prepared for creating X1_inclusive;
for i = 1:Nr
    if isnan(RAW{i,1})
        continue;
    end
    tmp1 = strsplit(RAW{i,1},',');
    tmp1 = strrep(tmp1,' ','');
    level1_variant{i,1} = tmp1{1};  %variant of this level;
    for j = 1:numel(tmp1)
        ind = strcmp(variant_name,tmp1{j});
        X1tmp(i,barcode(ind,:)==1) = 1;
    end
end
level1_variant = level1_variant(~cellfun(@isempty,level1_variant));
variant_level1_only = unique(level1_variant);
Nlevel1 = numel(variant_level1_only);


X1_inclusive = zeros(Nlevel1,Nsnp);  %X for level1-variant and its sublineages; from covspectrum;
X1_only = zeros(Nlevel1,Nsnp); %X for level1-variant only from covspectrum; should be a subset of X1tmp;
variant_level1_inclusive = cell(Nlevel1,1);
for i = 1:(Nlevel1)
    ind = find(strcmp(level1_variant,variant_level1_only{i}));
    tmp1 = X1tmp(ind,:);
    if numel(ind) == 1
        X1_inclusive(i,(tmp1)>0) = 1;
    else
        X1_inclusive(i,max(tmp1)>0) = 1;
    end
    ind2 = find(strcmp(variant_name,variant_level1_only{i}));
    X1_only(i,:) = barcode(ind2,:);
    ind3 = find(contains(RAW,variant_level1_only{i}));

    for j = 1:numel(ind3)
        tmp1 = strsplit(RAW{ind3(j),1},',');
        tmp1 = strrep(tmp1,' ','');
        if j == 1
            variant_level1_inclusive{i,1} = tmp1(:);
        else
            variant_level1_inclusive{i,1} = [variant_level1_inclusive{i,1};tmp1(:)];
        end
    end
    variant_level1_inclusive{i,1} = unique(variant_level1_inclusive{i,1});
end
end