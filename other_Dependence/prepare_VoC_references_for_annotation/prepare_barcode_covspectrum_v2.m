function prepare_barcode_covspectrum_v2(family_tree_file,cov_spectrum_path,covspectrum_thre)
%% prepare X for each variant; using covspectrum all mutations;
% use hierachical structure:
%   level1: include all mutations with proportion>0.9;
%   level2: include all mutations with proportion>0.9 & remove those in
%           level1-parent;
%   level3: include all mutations with proportion>0.9 & remove those in
%           level1 and level2-parent; 
%   level4: include all mutations with proportion>0.9 & remove those in
%           level1 and level2 and level3 -parent; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath([cov_spectrum_path,'\..\common_func\'])  %determinant_mutations_for_each_lineage.m;
                                        % determine_barcode_for_each_hierachical_level.m
if ~exist('covspectrum_thre','var') 
covspectrum_thre = 0.9;  %only true ones belong to this variant; more than 90% clinical samples should contain this mutation
end
%% get all Variant of interest from hierachical tree
[~,~,RAW] = xlsread([cov_spectrum_path,'\',family_tree_file]);
[Nr,Nc] = size(RAW);
count = 1;
for i = 1:Nr
    for j = 1:Nc
        if isnan(RAW{i,j})
            continue;
        end
        tmp1 = strsplit(RAW{i,j},',');
        tmp1 = strrep(tmp1,' ','');
        N1 = numel(tmp1);
        Variant_of_interest(count:count+N1-1) = tmp1;
        count = count + N1;
    end
end
Variant_of_interest = unique(Variant_of_interest);
N_VOI = numel(Variant_of_interest);

%% prepare X_covspectrum: include all mutations in covspectrum.org with all-range frequencies; not threshold
Determinant_mutation_covspectrum = cell(N_VOI,1);
for i = 1:N_VOI
    T1 = readcell([cov_spectrum_path,'\nt\',Variant_of_interest{i},'.csv']);
    Determinant_mutation_covspectrum{i,1} = T1(2:end,1);
    Determinant_mutation_covspectrum{i,2} = cell2mat(T1(2:end,2)); %proportion
end
%% save cov-spectrum
for i = 1:N_VOI
    if i==1
        SNP_name = Determinant_mutation_covspectrum{i,1};
    else
        SNP_name = [SNP_name;Determinant_mutation_covspectrum{i,1}];
    end
end
SNP_name_f = unique(SNP_name); Nf = numel(SNP_name_f);
X_covspectrum = zeros(N_VOI,Nf);
for i = 1:N_VOI
    [~,ind1,ind2] = intersect(SNP_name_f,Determinant_mutation_covspectrum{i,1});
    X_covspectrum(i,ind1) = Determinant_mutation_covspectrum{i,2}(ind2);
end

%% threshold with frequecies > covspectrum_thre: 0.9; to include only determinant mutations
ind_pass_thre = find(sum(X_covspectrum<covspectrum_thre)<numel(Variant_of_interest));
barcode = double(X_covspectrum(:,ind_pass_thre)>=covspectrum_thre);
Variant_of_interest = Variant_of_interest(:); Nvariant_in_barcode = numel(Variant_of_interest);
loc_barcode = SNP_name_f(ind_pass_thre);  Nsnp = numel(loc_barcode);

%% level1 variant are the first element before ',' in RAW(:,1); include
% level1,2,3 lineages in X1_inclusive; and mutations in level1 variants only in X1_only;
[X1_inclusive,level1_variant_inclusive,X1_only,level1_variant] = ...
    determine_barcode_for_each_hierachical_level(RAW(:,1),barcode,Variant_of_interest); 

%% level2: include level2 and level3 sublineages;
RAW2 = RAW(:,2);
RAW2 = RAW2(XZ_non_nan_cell(RAW2));
[X2_inclusive,level2_variant_inclusive,X2_only,level2_variant] = ...
    determine_barcode_for_each_hierachical_level(RAW2,barcode,Variant_of_interest);
X2_inclusive_remove_level1 = X2_inclusive;
X2_only_remove_level1 = X2_only;
% remove those in X1_only;
count = 1;
for i = 1:numel(level2_variant)
    for j = 1:numel(level1_variant_inclusive)
        ind = find(contains(level1_variant_inclusive{j},level2_variant{i}));
        if ~isempty(ind)
            X2_only_remove_level1(i,X1_only(j,:)==1) = 0;
            X2_inclusive_remove_level1(i,X1_only(j,:) == 1) = 0;
            count = count + 1;
            break;
        end
    end
end
% determine parent level
level2_variant_parent = cell(numel(level2_variant),1);
for i = 1:numel(level2_variant)
    for j = 1:numel(level1_variant)
        ind = find(contains(level1_variant_inclusive{j},level2_variant{i}));
        if ~isempty(ind)
            level2_variant_parent{i,1} = level1_variant{j};
            break;
        end
    end
end
%% level3:
RAW3 = RAW(:,3);
RAW3 = RAW3(XZ_non_nan_cell(RAW3));
[X3_inclusive,level3_variant_inclusive,X3_only,level3_variant] = ...
    determine_barcode_for_each_hierachical_level(RAW3,barcode,Variant_of_interest);
X3_inclusive_remove_level12 = X3_inclusive;
X3_only_remove_level12 = X3_only;
% remove those in X1_only and X2_only;
for i = 1:numel(level3_variant)
    for j = 1:numel(level1_variant_inclusive)
        ind = find(contains(level1_variant_inclusive{j},level3_variant{i}));
        if ~isempty(ind)
            X3_only_remove_level12(i,X1_only(j,:)==1) = 0;
            X3_inclusive_remove_level12(i,X1_only(j,:) == 1) = 0;
            break;
        end
    end
    for j = 1:numel(level2_variant_inclusive)
        ind = find(contains(level2_variant_inclusive{j},level3_variant{i}));
        if ~isempty(ind)
            X3_only_remove_level12(i,X2_only(j,:)==1) = 0;
            X3_inclusive_remove_level12(i,X2_only(j,:) == 1) = 0;
            break;
        end
    end
end
% determine parent level
level3_variant_parent = cell(numel(level3_variant),2);
for i = 1:numel(level3_variant)
    for j = 1:numel(level1_variant)
        ind = find(contains(level1_variant_inclusive{j},level3_variant{i}));
        if ~isempty(ind)
            level3_variant_parent{i,1} = level1_variant{j};
            break;
        end
    end
    for j = 1:numel(level2_variant)
        ind = find(contains(level2_variant_inclusive{j},level3_variant{i}));
        if ~isempty(ind)
            level3_variant_parent{i,2} = level2_variant{j};
            break;
        end
    end
end


%% level4:
RAW4 = RAW(:,4);
RAW4 = RAW4(XZ_non_nan_cell(RAW4));
[X4_inclusive,level4_variant_inclusive,X4_only,level4_variant] = ...
    determine_barcode_for_each_hierachical_level(RAW4,barcode,Variant_of_interest);
X4_inclusive_remove_level123 = X4_inclusive;
X4_only_remove_level123 = X4_only;
% remove those in X1_only X2_only and X3_only;
for i = 1:numel(level4_variant)
    for j = 1:numel(level1_variant_inclusive)
        ind = find(contains(level1_variant_inclusive{j},level4_variant{i}));
        if ~isempty(ind)
            X4_only_remove_level123(i,X1_only(j,:)==1) = 0;
            X4_inclusive_remove_level123(i,X1_only(j,:) == 1) = 0;
            break;
        end
    end
    for j = 1:numel(level2_variant_inclusive)
        ind = find(contains(level2_variant_inclusive{j},level4_variant{i}));
        if ~isempty(ind)
            X4_only_remove_level123(i,X2_only(j,:)==1) = 0;
            X4_inclusive_remove_level123(i,X2_only(j,:) == 1) = 0;
            break;
        end
    end
    for j = 1:numel(level3_variant_inclusive)
        ind = find(contains(level3_variant_inclusive{j},level4_variant{i}));
        if ~isempty(ind)
            X4_only_remove_level123(i,X3_only(j,:)==1) = 0;
            X4_inclusive_remove_level123(i,X3_only(j,:) == 1) = 0;
            break;
        end
    end
end
% determine parent level
level4_variant_parent = cell(numel(level4_variant),2);
for i = 1:numel(level4_variant)
    for j = 1:numel(level1_variant)
        ind = find(contains(level1_variant_inclusive{j},level4_variant{i}));
        if ~isempty(ind)
            level4_variant_parent{i,1} = level1_variant{j};
            break;
        end
    end
    for j = 1:numel(level2_variant)
        ind = find(contains(level2_variant_inclusive{j},level4_variant{i}));
        if ~isempty(ind)
            level4_variant_parent{i,2} = level2_variant{j};
            break;
        end
    end
    for j = 1:numel(level3_variant)
        ind = find(contains(level3_variant_inclusive{j},level4_variant{i}));
        if ~isempty(ind)
            level4_variant_parent{i,3} = level3_variant{j};
            break;
        end
    end
end

%%======================================================================================
barcode_f = zeros(Nvariant_in_barcode,Nsnp);
barcode_level = zeros(Nvariant_in_barcode,1);
barcode_parent = cell(Nvariant_in_barcode,3);
for i = 1:Nvariant_in_barcode
    ind1 = find(strcmp(level1_variant,Variant_of_interest{i}));
    if ~isempty(ind1)  %level1; use all determinant mutations > 0.9; but only itself; not inclusive;
        barcode_level(i,1) = 1;
        barcode_f(i,:) = X1_only(ind1,:);
        continue;
    end
    ind2 = find(strcmp(level2_variant,Variant_of_interest{i}));
    if ~isempty(ind2)  %level2; use all determinant mutations > 0.9 after removing those in level1;
        barcode_level(i,1) = 2;
        barcode_f(i,:) = X2_only_remove_level1(ind2,:);
        barcode_parent(i,1) = level2_variant_parent(ind2,1); 
        continue;
    end
    ind3 = find(strcmp(level3_variant,Variant_of_interest{i}));
    if ~isempty(ind3)  %level3; use all determinant mutations > 0.9 after removing those in level1 and level2;
        barcode_level(i,1) = 3;
        barcode_f(i,:) = X3_only_remove_level12(ind3,:);
        barcode_parent(i,1:2) = level3_variant_parent(ind3,:); 
        continue;
    end

    ind4 = find(strcmp(level4_variant,Variant_of_interest{i}));
    if ~isempty(ind4)  %level3; use all determinant mutations > 0.9 after removing those in level1 and level2;
        barcode_level(i,1) = 4;
        barcode_f(i,:) = X4_only_remove_level123(ind4,:);
        barcode_parent(i,:) = level4_variant_parent(ind4,:);
    else
        disp('error-no-this-variant')
    end
end
[variant_determinane_mutationsf] = determinant_mutations_for_each_lineage(barcode_f,Variant_of_interest,loc_barcode);

xlsout_barcode = cell(Nvariant_in_barcode,7);
xlsout_barcode(:,1) = variant_determinane_mutationsf(:,1);
for i = 1:Nvariant_in_barcode
    mutation_tmp = variant_determinane_mutationsf{i,2};
    xlsout_barcode{i,2} = numel(mutation_tmp);
    for j = 1:xlsout_barcode{i,2}
        if j == 1
            tmp1 = mutation_tmp{j};
        else
            tmp1 = [tmp1, ',', mutation_tmp{j}];
        end
    end
    xlsout_barcode{i,7} = tmp1;
    clear tmp1;
end
xlsout_barcode(:,3) = num2cell(barcode_level);
xlsout_barcode(:,4:6) = barcode_parent;
xlsoutf = [{'variant-name','number-of-mutations-in-barcode_f-XZ','variant-hierachical-level','level1-parent','level2-parent','level3-parent','mutations'};...
    xlsout_barcode];
xlswrite([cov_spectrum_path,'\covspectrum_variant_of_interest_hierachical.xlsx'],xlsoutf,[num2str(covspectrum_thre),'_v2']);
N_in_each_variant = sum(barcode_f,2);
xlsout_N = [Variant_of_interest,num2cell(N_in_each_variant)];

save([cov_spectrum_path,'\covspectrum_variant_of_interest_hierachical_',num2str(covspectrum_thre),'_v2.mat'],"variant_determinane_mutationsf",...
    "loc_barcode","barcode_f","barcode_level","barcode_parent","Variant_of_interest",...
    "level1_variant","level2_variant","level3_variant","level1_variant_inclusive","level2_variant_inclusive",'xlsout_N');
%% save dominant mutation list for each VOI
save_path = [cov_spectrum_path,'\covspectrum_dominant_mutations'];
mkdir_XZ(save_path);
for j = 1:N_VOI
    writecell(variant_determinane_mutationsf{j,2},[save_path,'\',variant_determinane_mutationsf{j,1},'_v2.txt']);
end
