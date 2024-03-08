function [Yf_recode_del,feature_name_recode_del,feature_remove] = recode_deletion...
    (genome_fasta_file,feature_name,Yf)  
% keep
    % Yf: Nsampe x Nloc(Nf);
    % feature_name: Nloc(Nf) x 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin < 3
        Yf = [];
    end
%     T1 = readcell([root_path,'\covid_genome\wuhan.fasta'],'FileType','text','Delimiter','\t','NumHeaderLines',0);
    T1 = readcell(genome_fasta_file,'FileType','text','Delimiter','\t','NumHeaderLines',0);
    Nf = numel(feature_name);
    Nrow = size(T1);
    for i = 2:Nrow
        if i==2
            ref_fasta = T1{i,1};
        else
            ref_fasta = [ref_fasta T1{i,1}];
        end
    end
    disp(numel(ref_fasta)); %should be 29903;
    
    %% recode deletion in iVar: RefLoc-Deleted --> Ref@deleted[Loc+1]- (Ref@deleted = RefLoc+1;
    count = 1;
    count1 = 1; %number of insertions removed
    for i = 1:Nf
        if ~(contains(feature_name{i},'-') || contains(feature_name{i},'+'))
            feature_name_recode_del{count,1} = feature_name{i};
            if ~isempty(Yf)
                if ndims(Yf)==2
                    Yf_recode_del(:,count) = Yf(:,i);
                elseif ndims(Yf)==3
                    Yf_recode_del(:,:,count) = Yf(:,:,i);
                end
            end
            count = count + 1;
        else
            % insertion or deletion;
            tmp_feature = feature_name{i};
            tmp_loc = regexp(tmp_feature,'\d*','match'); 
            tmp1 = strsplit(tmp_feature,tmp_loc);
            Ntmp_indel = numel(tmp1{2});

            tmp_loc_num = str2double(tmp_loc);
            tmp_ref = ref_fasta(tmp_loc_num);
            if ~strcmp(tmp1,tmp_ref)
                disp('error: ref allele not match 1');
            end
            if strcmp(tmp1{2}(1),'-')
                % deletion; ivar coded A123-G is deleting G at 124
                % position;
                for j = 2:Ntmp_indel  %first is '-';
                    current_ref =ref_fasta(tmp_loc_num+j-1);
                    if strcmp(current_ref,tmp1{2}(j))
                        feature_name_recode_del{count,1} = [current_ref,num2str(tmp_loc_num+j-1),'-'];
                        if ~isempty(Yf)
                            if ndims(Yf)==2
                                Yf_recode_del(:,count) = Yf(:,i);
                            elseif ndims(Yf)==3
                                Yf_recode_del(:,:,count) = Yf(:,:,i);
                            end
                        end
                        count = count + 1;
                    else
                        disp('error: ref allele not match 2');
                    end
                end
            else
                % insertion still not consider
                feature_remove{count1,1} = tmp_feature;
                count1 = count1 + 1;
%                 disp(['remove',tmp_feature]);
            end
        end
    end
    if ~exist("Yf_recode_del",'var')
        Yf_recode_del =[];
    end
end