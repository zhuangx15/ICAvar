function [determinane_mutationsf] = determinant_mutations_for_each_lineage(X,variant_vec,loc_barcode)
    Nvariant = numel(variant_vec);
    determinant_mutations = cell(Nvariant,1);
    for i = 1:Nvariant
        ind = find(X(i,:) == 1);
        determinant_mutations{i,1} = loc_barcode(ind);
    end
    determinane_mutationsf = [variant_vec(:) determinant_mutations];
end