% use the mask distance as well as the main branch distance
% collect the latent label 
% eavluate the effect of auto context
% includes the latent label 
% discard the kernel features and adopts the new admm features
% samples_idx(:,1) => sample image no
% samples_idx(:,2) => sample row
% samples_idx(:,3) => sample column
% samples_idx(:,4) => sample label (-1/+1)
function [weak_learners_admm] = train_admm_lat_3D(params,fn,samples_idx)


samples_no = size(samples_idx,1);

labels = samples_idx(:,5);

samples_idx = samples_idx(:,1:4);


current_response_admm = zeros(samples_no,1);

current_response_ac = current_response_admm;

current_response_ctx = current_response_admm;


[compute_wi,compute_ri,compute_loss,compute_indiv_loss,compute_2nd_deriv,mex_loss_type] = select_fncts(params,labels);

params.mex_loss_type = mex_loss_type; 


W = compute_wi(current_response_admm);

R = compute_ri(current_response_admm);

W_ac = W;

R_ac = R;

W_ctx = W;

R_ctx = R;


% assign the weight for each individual data point

% wgt_samp = weight_sample_tol_3D(fn.train.gts,params,samples_idx);


wgt_samp = ones(size(samples_idx,1),1);


% X = data.train.imgs.X(:,data.train.imgs.idxs);

[features_admm,ftrs_kernel,ftrs_win] = collect_admm_ftrs_3D(fn.train.imgs,samples_idx);

% [features_admm1,ftrs_kernel1,ftrs_win1] = collect_admm_ftrs1_3D(fn.train.imgs,samples_idx);

% features_admm = [features_admm,features_admm1];

% ftrs_kernel = [ftrs_kernel;ftrs_kernel1];

% ftrs_win = [ftrs_win;ftrs_win1];

for i_w = 1:params.wl_no
    
    t_wl = tic;
    
    fprintf('  Learning WL %d/%d\n',i_w,params.wl_no);
    
    % Indexes of the two training subparts
    T1_idx = sort(randperm(length(labels),params.T1_size),'ascend');
    T2_idx = sort(randperm(length(labels),params.T2_size),'ascend');
    [wr_idxs,wr_responses,wr_weights] = compute_wr(params,T1_idx,W,R,compute_indiv_loss,compute_2nd_deriv,labels,current_response_admm);

    
    W = W .* wgt_samp;
    
    W_ctx = W_ctx .* wgt_samp;
    
    W_ac = W_ac .* wgt_samp;
    
    
    
    %%%%%%%%%%%%%%%%% train the base line weak learners%%%%%%%%%%%%%%%%%
    
%     
%     [weak_learners(i_w).kernels,weak_learners(i_w).kernel_params,...
%         weak_learners(i_w).reg_tree] = train_admm_reg(X,params,samples_idx,T1_idx,T2_idx,W,R);
    
    
   % collect the existing kernel boost features 
    
%     [kernels_kb,kernel_params_kb,features_kb] = train_kernel_features(X,params,samples_idx,T1_idx,T2_idx,W,R);
    
    % collect the admm features
    
    [weak_learners_admm(i_w).kernels,weak_learners_admm(i_w).kernel_params,...
        weak_learners_admm(i_w).reg_tree,weak_learners_admm(i_w).admm_list] =...
        train_admm_reg_3D(params,features_admm,T2_idx,W,R,ftrs_kernel,ftrs_win);
    
    fprintf('    Performing prediction on the whole training set...\n');
    
    admm_list = weak_learners_admm(i_w).admm_list;
    
    t_pr = tic;
    cached_responses_admm = LDARegStumpPredict(weak_learners_admm(i_w).reg_tree,...
        single(features_admm(:,admm_list)));
    time_pr = toc(t_pr);
    fprintf('    Prediction finished, took %f seconds\n',time_pr);
    
    weak_learners_admm(i_w).alpha = search_alpha(current_response_admm,cached_responses_admm,labels,params);
    
    current_response_admm = current_response_admm + weak_learners_admm(i_w).alpha * cached_responses_admm;
    
    
    W = compute_wi(current_response_admm);
    
    R = compute_ri(current_response_admm);

    train_scores(i_w,:) = weak_learner_scores(current_response_admm,labels,wgt_samp,compute_loss);
   
    
    save([params.codename '_weak_learners_sav.mat'],'weak_learners_admm');
    
    save([params.codename '_train_scores_sav.mat'],'train_scores');
    
    
    
    %%%%%%%%%%%%%%%%%%%%  complete training the base line weak learners %%
    
    wl_time = toc(t_wl);
    
    fprintf('  Learning WL %d took %f seconds\n------------------------------------------------\n\n',i_w,wl_time);
    
end

end
