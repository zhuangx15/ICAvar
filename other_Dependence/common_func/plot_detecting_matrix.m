function [h3,xlsoutf] = plot_detecting_matrix(detect_matrix_ICA_corr,detect_matrix_ICA_corr_binary,date_vec,ind_date_display,...
    variant_name,cmap,matrix_type_string,ICA_criteria_corr)
if nargin < 7
    matrix_type_string = 'ICA';
end
% xlsoutf: first and last detect date by this method;
Nvariant_in_barcode = numel(variant_name);
if isempty(detect_matrix_ICA_corr_binary)
    h3 = figure(randi(5000));set(h3,'Position',[50,50,3200,600])
    imagesc(detect_matrix_ICA_corr(:,ind_date_display));
    if ~isempty(cmap)
        colormap(cmap); clim([-1,1]); 
    else
        clim([0,1]);
    end
    clb = colorbar; clb.Label.String = matrix_type_string;
    xticks(1:numel(ind_date_display));xticklabels(date_vec(ind_date_display));
    yticks(1:Nvariant_in_barcode);yticklabels(variant_name);
    title(['detected-by-',matrix_type_string]);
    box on; grid on;
    set(gca,'FontSize',16);
    set(gca,'FontWeight','bold');
    xlsoutf = [];
    
elseif isempty(detect_matrix_ICA_corr)
    h3 = figure(randi(5000));set(h3,'Position',[50,50,3200,600])
    imagesc(detect_matrix_ICA_corr_binary(:,ind_date_display));
    if ~isempty(cmap)
        colormap(cmap); clim([-1,1]); 
    else
        clim([0,1]);
    end

    clb = colorbar; clb.Label.String = matrix_type_string;
    xticks(1:numel(ind_date_display));xticklabels(date_vec(ind_date_display));
    yticks(1:Nvariant_in_barcode);yticklabels(variant_name);
    title(['detected-by-',matrix_type_string]);
    box on; grid on;
    set(gca,'FontSize',16);
    set(gca,'FontWeight','bold');
    [xlsout_detect_date_ICA] = determine_first_and_last_Detect_date(detect_matrix_ICA_corr_binary(:,ind_date_display),...
        date_vec(ind_date_display),variant_name);
    xlsoutf = [{'variant','first-ICA','last-ICA'};xlsout_detect_date_ICA(:,1:3)];
else
    h3 = figure(randi(5000));set(h3,'Position',[50,50,3200,600])
    subplot(121);
    imagesc(detect_matrix_ICA_corr(:,ind_date_display));
    if ~isempty(cmap)
        colormap(cmap); clim([-1,1]); 
    else
        clim([0,1]);
    end
    clb = colorbar; clb.Label.String = matrix_type_string;
    xticks(1:numel(ind_date_display));xticklabels(date_vec(ind_date_display));
    yticks(1:Nvariant_in_barcode);yticklabels(variant_name);
    title(['detected-by-',matrix_type_string]);
    box on; grid on;
    subplot(122);
    imagesc(detect_matrix_ICA_corr_binary(:,ind_date_display));
    if ~isempty(cmap)
        colormap(cmap); clim([-1,1]); 
    else
        clim([0,1]);
    end
    xticks(1:numel(ind_date_display));xticklabels(date_vec(ind_date_display));
    yticks(1:Nvariant_in_barcode);yticklabels(variant_name);
    title(ICA_criteria_corr);
    box on; grid on;
    set(gca,'FontSize',16);
    set(gca,'FontWeight','bold');
    [xlsout_detect_date_ICA] = determine_first_and_last_Detect_date(detect_matrix_ICA_corr_binary(:,ind_date_display),...
        date_vec(ind_date_display),variant_name);
    xlsoutf = [{'variant','first-ICA','last-ICA'};xlsout_detect_date_ICA(:,1:3)];
end
end


function [xlsout_detect_date] = determine_first_and_last_Detect_date(detect_matrix_ICA,date_vec,variant_name) 
%detect_matrix: Nvariant x Ndate; binarized matrix;

Nvariant_in_barcode = numel(variant_name);
xlsout_detect_date = cell(Nvariant_in_barcode,3);
for i = 1:Nvariant_in_barcode
    ind1 = find(~isnan(detect_matrix_ICA(i,:)) & detect_matrix_ICA(i,:)==1);
    if ~isempty(ind1)
        xlsout_detect_date{i,1} = date_vec{ind1(1)};
        xlsout_detect_date{i,2} = date_vec{ind1(end)};
    end
end
xlsout_detect_date = [variant_name(:) xlsout_detect_date];
end