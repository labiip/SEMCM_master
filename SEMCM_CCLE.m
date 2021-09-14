clear all;close all;clc;

%addpath('SEMCM-master');

DEBUG = 0;
CCLE = 0;

% import matrix for completion

Y_filepath = 'ccle.xlsx';
[Y, text, alldata] = xlsread(Y_filepath);
Y_dim = size(Y);

Known_index0 = find(isnan(Y)==0);
Missing_index0 = find(isnan(Y)==1);
missing_index_1 = randperm(length(Known_index0),500);
Known_index = setdiff(Known_index0,Known_index0(missing_index_1));
Missing_index = union(Missing_index0,Known_index0(missing_index_1));

% replace NaN with 0 in Y
Y_0 = Y;
Y_0(Missing_index) = 0;

Mask_mat = ones(size(Y));
Mask_mat(Missing_index) = 0;
Mask_mat = logical(Mask_mat);
    
Y_recover1 = zeros(Y_dim);
Y_recover2 = Y_recover1;  
result = zeros(30,9);
for cv_run = 1:30
    rand('seed',cv_run + 20000);
    K_fold = 10;
    indices = crossvalind('Kfold',Known_index,K_fold);

for fold_id = 1:10
    % cross valind for Y
    fprintf('Fold %d:\n',fold_id);
    test = (indices == fold_id); train = ~test;
    Y_tmp = Y_0;
    Y_tmp(Known_index(test)) = 0;
%     max_Y_tmp = max(max(Y_tmp));
%     mean_Y_tmp = (max(Y_tmp(Known_index(train)))-min(Y_tmp(Known_index(train))))/2;
%     Y_tmp = Y_tmp - mean_Y_tmp;
    % Y_tmp = Y_tmp / max_Y_tmp;
    Mask_test = Mask_mat;
    Mask_test(Known_index(test)) = 0;

    C_k = zeros(Y_dim(2));  % 491*491
    A_k = zeros(Y_dim(2));  % 491*491
    Y_k = zeros(Y_dim);     % 23*491
    X_k = zeros(Y_dim);     % 23*491
    E_k = zeros(Y_dim);     % 23*491
    Q1_k = zeros(Y_dim);    % 23*491
    Q2_k = zeros(Y_dim);    % 23*491
    Q3_k = zeros(Y_dim(2)); % 491*491
    obj_rec = [];
    max_iter = 1000;
    
    mu = 0.4;
    mu_max = 1e7;
    lambda = 8;
    rho = 1.01;
    alpha = 1.6;
    iter_id = 1;
    decrease_stat = 1;
    
    while iter_id < max_iter
        % update A_k
        if DEBUG >= 1
            disp('Updating A');
        end
        % A_k1 = inv(X_k.'*X_k + eye(Y_dim(2)))*(X_k.'*(X_k - E_k + Q1_k/mu)+C_k+Q3_k/mu);
        A_k1 = (X_k.'*X_k + eye(Y_dim(2)))\(X_k.'*(X_k - E_k + Q1_k/mu)+C_k+Q3_k/mu);
        
        % update C_k
        if DEBUG >= 1
            disp('Updating C');
        end
        tmp = A_k1-Q3_k/mu;
        C_k1 = max(abs(tmp)-1/mu,0).*sign(tmp);
        % C_k1 = C_k1 - diag(C_k1);
        C_k1(logical(eye(size(C_k1)))) = 0;
        
        % update X_k
        if DEBUG >= 1
            disp('Updating X');
        end
%         X_k1 = ((E_k - Q1_k/mu)*(eye(Y_dim(2))-A_k1).' + Y_k +Q2_k/mu)...
%                *inv((eye(Y_dim(2))-A_k1)*(eye(Y_dim(2))-A_k1).'+eye(Y_dim(2)));
        X_k1 = ((E_k - Q1_k/mu)*(eye(Y_dim(2))-A_k1).' + Y_k +Q2_k/mu)...
                /((eye(Y_dim(2))-A_k1)*(eye(Y_dim(2))-A_k1).'+eye(Y_dim(2)));
        X_k1(Mask_test) = Y_tmp(Mask_test);
        
        % update Y_k
        if DEBUG >= 1
            disp('Updating Y');
        end
        tmp = X_k1 - Q2_k/mu;
        [U,Sigma,V] = svd(tmp);
        tmp1 = diag(Sigma) - alpha/mu; 
        Sigma = diag(tmp1);
        Sigma = [Sigma; zeros(Y_dim(1)-Y_dim(2),Y_dim(2))];
        Y_k1 = U*Sigma*V.';
        
        % update E_k
        if DEBUG >= 1
            disp('Updating E');
        end
        % E_k1 = (2*lambda+mu)*(X_k1-X_k1*A_k1+Q1_k/mu)/mu;
        tmp = X_k1 - X_k1*A_k1 + Q1_k/mu;
        [U,Sigma,V] = svd(tmp);
        tmp1 = diag(Sigma) - lambda/mu; 
        Sigma = diag(tmp1);
        Sigma = [Sigma; zeros(Y_dim(1)-Y_dim(2),Y_dim(2))];
        E_k1 = U*Sigma*V.';
        
        % update Q_ks
        if DEBUG >= 1
            disp('Updating Qs');
        end
        Q1_k1 = Q1_k + mu*(X_k1-X_k1*A_k1-E_k1);
        Q2_k1 = Q2_k + mu*(Y_k1-X_k1);
        Q3_k1 = Q3_k + mu*(C_k1-A_k1);
        
        % update mu
        if DEBUG >= 1
            disp('Updating mu');
        end
        mu = min(rho*mu, mu_max);
        
        % check convergence
        obj = norm(C_k1,1)+lambda*norm(E_k1,1)+trace(Q1_k1'*(X_k1-X_k1*A_k1-E_k1))...
            + trace(Q2_k1'*(Y_k1-X_k1))+trace(Q3_k1*(C_k1-A_k1))...
            + mu/2*norm((X_k1-X_k1*A_k1-E_k1),2)+norm((Y_k1-X_k1),2);
        obj_rec = [obj_rec obj];
        
        if iter_id >= 2
            if obj_rec(end) < obj_rec(end - 1)
                if decrease_stat > 200 && abs(obj_rec(end) - obj_rec(end - 1))/obj_rec(end - 1) < 1e-5
                    break;
                else
                    decrease_stat = decrease_stat + 1;
                end
            else
                decrease_stat = 0;
            end
        end
        
        
        iter_id = iter_id + 1;
        C_k = C_k1;      % 491*491
        A_k = A_k1;      % 491*491
        Y_k = Y_k1;      % 23*491
        X_k = X_k1;      % 23*491
        E_k = E_k1;      % 23*491
        Q1_k = Q1_k1;    % 23*491
        Q2_k = Q2_k1;    % 23*491
        Q3_k = Q3_k1;    % 491*491
       
        
    end
  
   % Y_k = Y_k * max_Y_tmp;
%    Y_k = Y_k + mean_Y_tmp;
   Y_recover1(Known_index(test)) = Y_k(Known_index(test));
   Y_recover2(Known_index(test)) = X_k(Known_index(test));
end

figure;
subplot(2,2,1);imagesc(Y_0);colorbar;
subplot(2,2,2);imagesc(Y_recover1);colorbar;
subplot(2,2,3);scatter(Y_0(Known_index),Y_recover1(Known_index),'.');%corr(Y_0(Known_index),Y_recover1(Known_index))
subplot(2,2,4);plot(1:length(obj_rec),obj_rec);
% scatter(Y_0(Known_index),Y_recover2(Known_index),'.');%corr(Y_0(Known_index),Y_recover2(Known_index))
global_pcc = corr(Y_0(Known_index),Y_recover1(Known_index));
num = Y_0';
numpred = Y_recover1';

drugwisecorr = NaN(size(num,1),1);
drugwise_qt = NaN(size(num,1),1);
drugwiseerr = NaN(size(num,1),1);
drugwiseerr_qt = NaN(size(num,1),1);
drugwiserepn = NaN(size(num,1),1);
for d = 1:size(num,1)
    curtemp1 = num(d,:);
    y1 = prctile(curtemp1,75);
    xia1 = find(curtemp1 >= y1);
    y2 = prctile(curtemp1,25);
    xia2 = find(curtemp1 <= y2);
    xia = [xia1,xia2];
    drugwise_qt(d) = corr(curtemp1(xia)',numpred(d,xia)');
    drugwiseerr_qt(d) = sqrt(sum((curtemp1(xia)-numpred(d,xia)).^2)/sum(~isnan(curtemp1(xia))));
    curtemp2 = numpred(d,:);
    curtemp2(isnan(curtemp1)) = [];
    curtemp1(isnan(curtemp1)) = [];
    drugwiserepn(d) = length(curtemp1);  
    drugwisecorr(d) = corr(curtemp1',curtemp2');
    drugwiseerr(d) = sqrt(sum((curtemp1-curtemp2).^2)/sum(~isnan(curtemp1)));
end
ave_pcc = mean(drugwisecorr);
std_pcc = std(drugwisecorr);
ave_err = mean(drugwiseerr);
std_err = std(drugwiseerr);
ave_pcc_sr = mean(drugwise_qt);
std_pcc_sr = std(drugwise_qt);
ave_err_sr = mean(drugwiseerr_qt);
std_err_sr = std(drugwiseerr_qt);
result(cv_run, :) = [global_pcc ave_pcc ave_pcc_sr ave_err ave_err_sr std_pcc std_pcc_sr std_err std_err_sr];
end
