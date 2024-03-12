function [R_square_model,Bfull,regressor_T,regressor_P,regressor_D,...
    con_con_map,con_t_map,con_p_map,con_d_map,rho_map,var_con_con_map,...
    RES,dfe,Xfull,SSE,SSR,SST,F_model,p_model,ind_empty_col_in_X] = XZ_SPM_univariate_analysis...
    (X,Y,con,flag_add_constant,nuisance_in_model_flag_f)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% !!!multiple linear regression!!!
% !!!Xfull_inv = inv(Xfull'*Xfull); %pxp; diagonal: px1;  Here should us inv instead of pinv
% !!!   results using inv is consistent with glmfit; using pinv is not;

% input preprocessed data; if need z-score then z-scored already
% data: t*q; Number of voxel; multiple voxels and multiple regreesors in linear regression
% X: t*p;  Number of regressor
% con: N_con * p
%%% same Bfull, regressor_T, regressor_P, regressor_D results as glmfit !!tested!!
%%% add Confidence interval for each Beta; 09/08/2020

%%% remove all zero columns in DesignMatrix, and computes contrast only
%%% when sum(constrast)==0 after removal; 05252023
%%% add nuisance_in_model_flag_f to indicate if regressor in X are nuisance
%%% regressors, and remove them when computing correlation maps betwee (Y,
%%% Y_est): 05252023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 5 
    nuisance_in_model_flag_f = zeros(1,size(X,2));
end
if nargin  < 4
    flag_add_constant = 0;
end
if nargin < 3
    con = [];
end
N_con = size(con,1);

[tdim,qmaxM] = size(Y);  %txq;
if flag_add_constant == 1
    Xfull = [ones(tdim,1),X];
    con = [zeros(N_con,1),con];
    nuisance_in_model_flag_f = [0,nuisance_in_model_flag_f];
else
    Xfull = X;
    Y = zscore(Y);
end
%%===============change 05252023: remove empty column in designMatrix=======
ind_empty_col_in_X = (sum(Xfull==0) == size(Xfull,1));
Xfull = Xfull(:,~ind_empty_col_in_X);
if ~isempty(con)
    con = con(:,~ind_empty_col_in_X);
end
nuisance_in_model_flag_f = nuisance_in_model_flag_f(:,~ind_empty_col_in_X);
%%=========================================================================
dofR = rank(Xfull)-1;  
Np = size(Xfull,2);
dfe = tdim - dofR - 1;  % Xfull is always final X for matrix inverse; and dfe is always tdim-rank(X_full);
Bfull = pinv(Xfull)*Y; %p*q 
Yest = Xfull*Bfull; %t*q
RES = Y-Yest;  %mean square error at each voxel;  %2357*236657
SSE = sum(RES.^2);%1xq
SST = sum((Y-(ones(tdim,1)*mean(Y))).^2); %1xq
SSR = SST - SSE;
F_model = (SSR./dofR)./(SSE./dfe);
p_model = Fval_2_pval(F_model,dofR, dfe);
R_square_model = 1-SSE./SST; %1xq
Xfull_inv = inv(Xfull'*Xfull); %pxp; diagonal: px1;  Here should us inv instead of pinv
varRES = SSE./dfe; %1xq  % residual ~ N(0,varRes); 
                         % each beta follows normal distribution Gaussian(true_beta,varRES*diag(inv(Xfull'X)))
                         % CORRECT checked 09/08/2020
var_Bfull = diag(Xfull_inv)*varRES; %pxq; 
regressor_T = Bfull ./sqrt(var_Bfull);  %significance of every regressor 
regressor_P = tval_2_pval(regressor_T,dfe,'both');
regressor_D = Bfull./sqrt((ones(Np,1)*varRES));  %regressor effect size
% regressor_E = regressor_T./(sqrt(regressor_T.^2 + dfe)); %09082020: regressor correlation after controlling cov;

%% change 05252023: remove motion in computing rho-map;
Y_est_model = Xfull(:,nuisance_in_model_flag_f==0) * Bfull(nuisance_in_model_flag_f==0,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho_map = zeros(qmaxM,1);
for i = 1:qmaxM
    rho_map(i) = corr(Y(:,i),Y_est_model(:,i));
end

% %% regressor
% con_vec_for_regressor = eye(Np,Np);
% regressor_T = zeros(Np,qmaxM);
% regressor_P = zeros(Np,qmaxM);
% regressor_D = zeros(Np,qmaxM);
% for i = 1:Np
%     C = con_vec_for_regressor(i,:);
%     se = sqrt(varRES.*(C*pinv(Xfull'*Xfull)*C')); %1xq
%     con_con_map = C * (Bfull);   %1*234456
%     regressor_T(i,:) = con_con_map./(se+eps);
%     regressor_P(i,:) = 2 * tcdf(-abs(regressor_T(i,:)), dfe);
%     regressor_D(i,:) = con_con_map./sqrt(varRES);
% end

%% contrast
con_t_map = zeros(N_con,qmaxM);
con_d_map = zeros(N_con,qmaxM);
con_con_map = zeros(N_con,qmaxM);
var_con_con_map = zeros(N_con,qmaxM);
con_p_map = zeros(N_con,qmaxM);
if ~isempty(con)
    %% for contrast;
    for i = 1:N_con
        C = con(i,:); %1xp
%         if abs(sum(C(:)))<1e-10
            se = sqrt(varRES.*(C*pinv(Xfull'*Xfull)*C')); %1xq
            con_con_map(i,:) = C * (Bfull);   %1*234456
            var_con_con_map(i,:) = se.^2;
            con_t_map(i,:) = con_con_map(i,:)./(se+eps);
            con_p_map(i,:) = 2 * tcdf(-abs(con_t_map(i,:)), dfe); 
            con_d_map(i,:) = con_con_map(i,:)./sqrt(varRES);
%         end
    end
end
end