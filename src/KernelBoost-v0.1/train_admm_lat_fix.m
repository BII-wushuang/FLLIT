% use the mask distance as well as the main branch distance
% collect the latent label 
% eavluate the effect of auto context
% includes the latent label 
% discard the kernel features and adopts the new admm features
% samples_idx(:,1) => sample image no
% samples_idx(:,2) => sample row
% samples_idx(:,3) => sample column
% samples_idx(:,4) => sample label (-1/+1)
function [weak_learners,weak_learners_admm] = train_admm_lat_fix(params,data,samples_idx)


samples_no = size(samples_idx,1);

labels = samples_idx(:,4);

samples_idx = samples_idx(:,1:3);


current_response = zeros(samples_no,1);

current_response_ac = current_response;

current_response_ctx = current_response;


[compute_wi,compute_ri,compute_loss,compute_indiv_loss,compute_2nd_deriv,mex_loss_type] = select_fncts(params,labels);

params.mex_loss_type = mex_loss_type; 


W = compute_wi(current_response);

R = compute_ri(current_response);

W_ac = W;

R_ac = R;

W_ctx = W;

R_ctx = R;


% assign the weight for each individual data point

wgt_samp = weight_sample_tol(data.train.gts,params,samples_idx);


X = data.train.imgs.X(:,data.train.imgs.idxs);

features_admm = collect_admm_ftrs(X,samples_idx);

features_admm1 = collect_admm_ftrs1(X,samples_idx);

features_admm = [features_admm,features_admm1];

T1_idx = sort(randperm(length(labels),params.T1_size),'ascend');


i_ch = 1;

ch = 'imgs';

sub_ch_no = data.train.imgs.sub_ch_no;



fprintf('      Learning filters on the sub-channels\n');
for i_s = 1:sub_ch_no
    
    t_sch = tic;
    fprintf('        Learning on subchannel %d/%d of channel %s\n',i_s,sub_ch_no,ch);
    
    [kernels{i_ch}{i_s},kernel_params{i_ch}{i_s}] = mexMultipleSmoothRegression...
        (params,params.(ch),X(:,i_s),1 : sub_ch_no,samples_idx(T1_idx,1:3),R(T1_idx),W(T1_idx),i_ch,i_s,ch);
    
    sch_time = toc(t_sch);
    fprintf('        Completed, learned %d filters in %f seconds\n',length(kernels{i_ch}{i_s}),sch_time);
    
end



for i_s = 1 : sub_ch_no
    
    if(isfield(params,'n_kb'))
        
        kernels{i_ch}{i_s} = kernels{i_ch}{i_s}(1:params.n_kb);
        
        kernel_params{i_ch}{i_s} = kernel_params{i_ch}{i_s}(1:params.n_kb);
        
    end
    
end


for i_w = 1:params.wl_no
    
    t_wl = tic;
    
    fprintf('  Learning WL %d/%d\n',i_w,params.wl_no);
    
    % Indexes of the two training subparts
    T1_idx = sort(randperm(length(labels),params.T1_size),'ascend');
    T2_idx = sort(randperm(length(labels),params.T2_size),'ascend');
    [wr_idxs,wr_responses,wr_weights] = compute_wr(params,T1_idx,W,R,compute_indiv_loss,compute_2nd_deriv,labels,current_response);

    
    W = W .* wgt_samp;
    
    W_ctx = W_ctx .* wgt_samp;
    
    W_ac = W_ac .* wgt_samp;
    
    
    
    [weak_learners_admm(i_w).kernels,weak_learners_admm(i_w).kernel_params,...
        weak_learners_admm(i_w).reg_tree, weak_learners_admm(i_w).ctx_list]...
        = train_kernel_gb_ctx(X,kernels,kernel_params,params,features_admm,...
        samples_idx,W_ctx,R_ctx);
    
    cached_responses_ctx = evaluate_weak_learners_ctx(X,params,features_admm,samples_idx,weak_learners_admm(i_w));
   
    weak_learners_admm(i_w).alpha = search_alpha(current_response_ctx,cached_responses_ctx,labels,params);
    
    current_response_ctx = current_response_ctx + weak_learners_admm(i_w).alpha * cached_responses_ctx;    
    
    W_ctx = compute_wi(current_response_ctx);
    
    R_ctx = compute_ri(current_response_ctx);
    
    train_scores_ctx(i_w,:) = weak_learner_scores(current_response_ctx,labels,wgt_samp,compute_loss);    
    
    [weak_learners(i_w).kernels,weak_learners(i_w).kernel_params,...
        weak_learners(i_w).reg_tree, weak_learners(i_w).ctx_list]...
        = train_kernel_gb_ctx(X,kernels,kernel_params,params,[],...
        samples_idx,W,R);
    
    cached_responses = evaluate_weak_learners_ctx(X,params,[],samples_idx,weak_learners(i_w));
   
    weak_learners(i_w).alpha = search_alpha(current_response,cached_responses,labels,params);
    
    current_response = current_response + weak_learners(i_w).alpha * cached_responses;    
    
    W = compute_wi(current_response);
    
    R = compute_ri(current_response);
    
    train_scores(i_w,:) = weak_learner_scores(current_response,labels,wgt_samp,compute_loss);    
    
    save([params.codename '_weak_learners_sav.mat'],'weak_learners_admm','weak_learners');
    
    save([params.codename '_train_scores_sav.mat'],'train_scores_ctx','train_scores');
    
    %%%%%%%%%%%%%%%%%%%%  complete training the base line weak learners %%
    
    wl_time = toc(t_wl);
    
    fprintf('  Learning WL %d took %f seconds\n------------------------------------------------\n\n',i_w,wl_time);
    
end

end
