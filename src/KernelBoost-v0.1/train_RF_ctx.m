% use the mask distance as well as the main branch distance
% evaluate the effect of auto context
% takes the random forest framework

% samples_idx(:,1) => sample image no
% samples_idx(:,2) => sample row
% samples_idx(:,3) => sample column
% samples_idx(:,4) => sample label (-1/+1)
function [weak_learners,weak_learners_ctx] = train_RF_ctx(params,data,samples_idx)

% Train a KernelBoost classifier on the given samples


% specfically reserved for testing the effect of unchanged weight when
% extracting the features.


samples_no = size(samples_idx,1);

labels = samples_idx(:,4);

samples_idx = samples_idx(:,1:3);


current_response = zeros(samples_no,1);

current_response_ctx = current_response;


[compute_wi,compute_ri,compute_loss,compute_indiv_loss,compute_2nd_deriv,mex_loss_type] = select_fncts(params,labels);

params.mex_loss_type = mex_loss_type; 


W = compute_wi(current_response);

R = compute_ri(current_response);


W_ctx = W;

R_ctx = R;


% assign the weight for each individual data point

wgt_samp = weight_sample_tol(data.train.gts,params,samples_idx);


for i_w = 1:params.wl_no
    
    t_wl = tic;
    
    fprintf('  Learning WL %d/%d\n',i_w,params.wl_no);
    
    % Indexes of the two training subparts
    T1_idx = sort(randperm(length(labels),params.T1_size),'ascend');
    T2_idx = sort(randperm(length(labels),params.T2_size),'ascend');
    [wr_idxs,wr_responses,wr_weights] = compute_wr(params,T1_idx,W,R,compute_indiv_loss,compute_2nd_deriv,labels,current_response);

    
    W = W .* wgt_samp;
    
    
    % load the training images
    
    X = data.train.imgs.X(:,data.train.imgs.idxs);
    
    
    %%%%%%%%%%%%%%%%% train the base line weak learners%%%%%%%%%%%%%%%%%
    
    
    [weak_learners(i_w).kernels,weak_learners(i_w).kernel_params,...
        weak_learners(i_w).reg_tree] = train_kernel_boost(X,params,samples_idx,T1_idx,T2_idx,W,R);
    
    
    cached_responses = evaluate_weak_learners(X,params,samples_idx,weak_learners(i_w));
    
    weak_learners(i_w).alpha = search_alpha(current_response,cached_responses,labels,params);
    
    current_response = current_response + weak_learners(i_w).alpha * cached_responses;
    
    W = compute_wi(current_response);
    
    R = compute_ri(current_response);

    train_scores(i_w,:) = weak_learner_scores(current_response,labels,wgt_samp,compute_loss);
   
    
    %%%%%%%%%%%%%%%%%%%%  complete training the base line weak learners %%
    
end


% collect the global features as well as auto context features, be aware
% of the global features, they are only prepared for the next round, so
% this part should be behind the context and auto context learning
ctx_glb_all = collect_global_ctx(X,data.train.masks,...
    params,weak_learners(params.wl_no),samples_idx);


W_ctx = W_ctx .* wgt_samp;

for i_w = 1 : params.wl_no_ctx
    
    fprintf('  Learning WL %d/%d\n',i_w,params.wl_no);
    
    
    T1_idx = sort(randperm(length(labels),params.T1_size),'ascend');
    
    T2_idx = sort(randperm(length(labels),params.T2_size),'ascend');
    
    %%%%%%%%%%%% start training the ctx augmented weak learners%%%%%%%%%
    
    %  load the global context features
    
    
    [weak_learners_ctx(i_w).kernels,weak_learners_ctx(i_w).kernel_params,...
        weak_learners_ctx(i_w).reg_tree, weak_learners_ctx(i_w).ctx_list]...
        = train_kernel_rf_ctx(X,params,ctx_glb_all,...
        samples_idx,T1_idx,T2_idx,W_ctx,R_ctx);
    
    
    cached_responses_ctx = evaluate_weak_learners_ctx(X,params,ctx_glb_all,samples_idx,weak_learners_ctx(i_w));
   
    
    current_response_ctx = current_response_ctx + cached_responses_ctx;
    
    train_scores_ctx(i_w,:) = weak_learner_scores(current_response_ctx / i_w,labels,wgt_samp,compute_loss);
    
    
    save([params.codename '_weak_learners_sav.mat'],'weak_learners_ctx','weak_learners');
    
    save([params.codename '_train_scores_sav.mat'],'train_scores','train_scores_ctx');
    
    
    %%%% complete training the global context feature weak learners
    
    wl_time = toc(t_wl);
    
    fprintf('  Learning WL %d took %f seconds\n------------------------------------------------\n\n',i_w,wl_time);
    
end

end
