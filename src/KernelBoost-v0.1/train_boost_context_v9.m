% use the mask distance as well as the main branch distance
% samples_idx(:,1) => sample image no
% samples_idx(:,2) => sample row
% samples_idx(:,3) => sample column
% samples_idx(:,4) => sample label (-1/+1)
function [weak_learners,weak_learners_ctx1] = train_boost_context_v9(params,data,samples_idx)

% Train a KernelBoost classifier on the given samples
%
%  authors: Carlos Becker, Roberto Rigamonti, CVLab EPFL
%  e-mail: name <dot> surname <at> epfl <dot> ch
%  web: http://cvlab.epfl.ch/
%  date: February 2014

samples_no = size(samples_idx,1);

labels = samples_idx(:,4);
samples_idx = samples_idx(:,1:3);
current_response = zeros(samples_no,1);

[compute_wi,compute_ri,compute_loss,compute_indiv_loss,compute_2nd_deriv,mex_loss_type] = select_fncts(params,labels);

W = compute_wi(current_response);
R = compute_ri(current_response);

current_response_ctx1 = current_response;


W_ctx1 = W;

R_ctx1 = R;


n_km = 30;


% assign the weight for each individual data point

wgt_samp = zeros(size(samples_idx,1),1);

tol_z_wgt = params.tol_z_wgt;

fig_thres = params.fig_thres;


for i_img = 1 : length(data.train.gts.X)
    
    gt_img = data.train.gts.X{i_img};
    
    dist_gt = bwdist(gt_img);
    
    idx_img = (samples_idx(:,1) == i_img);
    
    samp_idx_img_2D = samples_idx(idx_img,2:3);
    
    samp_idx_img_1D = sub2ind(size(gt_img),samp_idx_img_2D(:,1),samp_idx_img_2D(:,2));
    
    wgt_img = dist_gt(samp_idx_img_1D);
    
    wgt_img(wgt_img > tol_z_wgt) = tol_z_wgt;

    wgt_img = wgt_img/ tol_z_wgt;
    
    wgt_img(dist_gt(samp_idx_img_1D) < 0.3) = 1;
    
    wgt_samp(idx_img) = wgt_img;
    
end

W = W .* wgt_samp;

W_ctx1 = W_ctx1 .* wgt_samp;

% W_ctx2 = W_ctx2 .* wgt_samp;


for i_w = 1:params.wl_no
    t_wl = tic;
    fprintf('  Learning WL %d/%d\n',i_w,params.wl_no);
    
    % Indexes of the two training subparts
    T1_idx = sort(randperm(length(labels),params.T1_size),'ascend');
    T2_idx = sort(randperm(length(labels),params.T2_size),'ascend');
    [wr_idxs,wr_responses,wr_weights] = compute_wr(params,T1_idx,W,R,compute_indiv_loss,compute_2nd_deriv,labels,current_response);    
    s_T1 = samples_idx(wr_idxs,1:3);
    s_T2 = samples_idx(T2_idx,1:3);
    
    features = cell(params.ch_no,1);
    kernels = cell(params.ch_no,1);
    kernel_params = cell(params.ch_no,1);
    
    for i_ch = 1:params.ch_no
        ch = params.ch_list{i_ch};
        fprintf('    Learning channel %s (%d/%d)\n',ch,i_ch,params.ch_no);
        X = data.train.(ch).X(:,data.train.(ch).idxs);
        X_idxs = data.train.(ch).idxs;
        
        sub_ch_no = data.train.(ch).sub_ch_no;
        features{i_ch} = cell(sub_ch_no,1);
        kernels{i_ch} = cell(sub_ch_no,1);
        kernel_params{i_ch} = cell(sub_ch_no,1);
        
        % Learn the filters
        fprintf('      Learning filters on the sub-channels\n');
        for i_s = 1:sub_ch_no
            t_sch = tic;
            fprintf('        Learning on subchannel %d/%d of channel %s\n',i_s,sub_ch_no,ch);
            [kernels{i_ch}{i_s},kernel_params{i_ch}{i_s}] = mexMultipleSmoothRegression(params,params.(ch),X(:,i_s),X_idxs,s_T1,wr_responses,wr_weights,i_ch,i_s,ch);
            sch_time = toc(t_sch);
            fprintf('        Completed, learned %d filters in %f seconds\n',length(kernels{i_ch}{i_s}),sch_time);
            
            t_ev = tic;
            fprintf('        Evaluating the filters learned on the subchannel\n');
            features{i_ch}{i_s} = mexEvaluateKernels(X(:,i_s),s_T2(:,1:3),params.sample_size,kernels{i_ch}{i_s},kernel_params{i_ch}{i_s});
            ev_time = toc(t_ev);
            fprintf('        Evaluation completed in %f seconds\n',ev_time);
        end
    end
    
    fprintf('    Merging features and kernels...\n');
    [kernels,kernel_params,features] = merge_features_kernels(kernels,kernel_params,features);
    fprintf('    Done!\n');
    

    
    
    
    fprintf('    Training regression tree on learned features...\n');
    t_tr = tic;
    reg_tree = LDARegStumpTrain(single(features),R(T2_idx),W(T2_idx)/sum(W(T2_idx)),uint32(params.tree_depth));
    time_tr = toc(t_tr);
    fprintf('    Done! (took %f seconds)\n',time_tr);
    
    
    if(i_w > 1)
        
        % from this step, the context featuers are added for individual image.
        
        nf = size(features,1);
        
        ncf_g = size(ctx_ftrs1,2);
        
        ctx_glb_all = zeros(size(samples_idx,1),ncf_g);
        
        
        for i_img = 1 : size(data.train.(ch).X,1)
            
            ctx_glb_all(samples_idx(:,1) == i_img,:) = ctx_ftrs1_w{i_w - 1,i_img};
            
            ctx_ftrs1_w{i_w - 1,i_img} = [];
            
        end
        
        ctx_glb = ctx_glb_all(T2_idx,:);

        
        fprintf('    Training regression tree on learned features combined with context 1...\n');
        t_tr = tic;
        
        if(isfield(params,'ctx1_tree_depth'))
              
            tree_d_g = params.ctx1_tree_depth;
              
        else
            
            tree_d_g = params.tree_depth;
            
        end
        
        reg_tree_g = LDARegStumpTrain(single([features,...
            ctx_glb]),R_ctx1(T2_idx),W_ctx1(T2_idx)...
            /sum(W_ctx1(T2_idx)),uint32(tree_d_g));
        time_tr = toc(t_tr);
        fprintf('    Done! (took %f seconds)\n',time_tr);

        
    end
    
    
    fprintf('    Removing useless kernels...\n');
    [weak_learners(i_w).kernels,weak_learners(i_w).kernel_params,weak_learners(i_w).reg_tree]...
        = remove_useless_filters(reg_tree,kernels,kernel_params);
    
    if(i_w > 1)
        
        [weak_learners_ctx1(i_w).kernels,weak_learners_ctx1(i_w).kernel_params,weak_learners_ctx1(i_w).reg_tree,...
            weak_learners_ctx1(i_w).ctx_list]...
            = remove_useless_filters_ctx(reg_tree_g,kernels,kernel_params);
        
        save([params.codename '_weak_learners_sav.mat'],'weak_learners_ctx1','weak_learners');
            
    end
    
    
    
    t_ev = tic;
    fprintf('    Evaluating the learned kernels on the whole training set...\n');
    features = zeros(length(labels),length(weak_learners(i_w).kernels));
    for i_ch = 1:params.ch_no
        ch = params.ch_list{i_ch};
        sub_ch_no = data.train.(ch).sub_ch_no;
        
        X = data.train.(ch).X(:,data.train.(ch).idxs);
        for i_s = 1:sub_ch_no
            idxs = find(cellfun(@(x)(x.ch_no==i_ch && x.sub_ch_no==i_s),weak_learners(i_w).kernel_params));
            if (~isempty(idxs))
                features(:,idxs) = mexEvaluateKernels(X(:,i_s),samples_idx(:,1:3),params.sample_size,weak_learners(i_w).kernels(idxs),weak_learners(i_w).kernel_params(idxs));
            end
        end
    end
    ev_time = toc(t_ev);
    fprintf('      Evaluation completed in %f seconds\n',ev_time);
    
    fprintf('    Performing prediction on the whole training set...\n');
    
    t_pr = tic;
    cached_responses = LDARegStumpPredict(weak_learners(i_w).reg_tree,single(features));
    time_pr = toc(t_pr);
    fprintf('    Prediction finished, took %f seconds\n',time_pr);
    
    clear features;
    
    if(i_w > 1)
        
        features_ctx1 = ctx_glb_all(:,weak_learners_ctx1(i_w).ctx_list);
        
        t_ev = tic;
        fprintf('    Evaluating the learned kernels on the whole training set...\n');
        features = zeros(length(labels),length(weak_learners_ctx1(i_w).kernels));
        for i_ch = 1:params.ch_no
            ch = params.ch_list{i_ch};
            sub_ch_no = data.train.(ch).sub_ch_no;
            
            X = data.train.(ch).X(:,data.train.(ch).idxs);
            for i_s = 1:sub_ch_no
                idxs = find(cellfun(@(x)(x.ch_no==i_ch && x.sub_ch_no==i_s),weak_learners_ctx1(i_w).kernel_params));
                if (~isempty(idxs))
                    features(:,idxs) = mexEvaluateKernels(X(:,i_s),...
                        samples_idx(:,1:3),params.sample_size,...
                        weak_learners_ctx1(i_w).kernels(idxs),weak_learners_ctx1(i_w).kernel_params(idxs));
                end
            end
        end
        ev_time = toc(t_ev);
        fprintf('      Evaluation completed in %f seconds\n',ev_time);
        
        fprintf('    Performing ctx1 prediction on the whole training set...\n');
        
        t_pr = tic;
        cached_responses_ctx1 = LDARegStumpPredict(...
            weak_learners_ctx1(i_w).reg_tree,single([features, features_ctx1]));
        time_pr = toc(t_pr);
        fprintf('    Prediction finished, took %f seconds\n',time_pr);
        
        clear features;
    
        clear features_ctx1 ctx_glb_all;
        
    end
    
    %clear features;
    
    fprintf('    Finding alpha through line search...\n');
    t_alp = tic;
    alpha = mexLineSearch(current_response,cached_responses,labels,mex_loss_type);
    time_alp = toc(t_alp);
    fprintf('    Good alpha found (alpha=%f), took %f seconds\n',alpha,time_alp);
    
    if(i_w > 1)
        
        fprintf('    Finding alpha 1 through line search...\n');
        t_alp = tic;
        alpha_ctx1 = mexLineSearch(current_response_ctx1,cached_responses_ctx1,labels,mex_loss_type);
        time_alp = toc(t_alp);
        fprintf('    Good alpha found (alpha=%f), took %f seconds\n',alpha_ctx1,time_alp);
        
        
    end
    
    
    alpha = alpha * params.shrinkage_factor;
    
    current_response = current_response + alpha*cached_responses;
    
    W = compute_wi(current_response);
    R = compute_ri(current_response);
    
    W = W .* wgt_samp;
    
    weak_learners(i_w).alpha = alpha;
    
    
    if(i_w > 1)
        
        alpha_ctx1 = alpha_ctx1 * params.shrinkage_factor;
        
        current_response_ctx1 = current_response_ctx1 + alpha_ctx1*cached_responses_ctx1;
        
        W_ctx1 = compute_wi(current_response_ctx1);
        
        R_ctx1 = compute_ri(current_response_ctx1);
        
        W_ctx1 = W_ctx1 .* wgt_samp;
        
        weak_learners_ctx1(i_w).alpha = alpha_ctx1;
       
        save([params.codename '_weak_learners_sav.mat'],'weak_learners','weak_learners_ctx1');
        
        
    end
    
    MR = sum((current_response>0)~=(labels>0))/length(labels);
    fprintf('    Misclassif rate: %.2f | Loss: %f\n',100*MR,compute_loss(current_response));
    
    train_scores(i_w,1) = 100*MR;
    train_scores(i_w,2) = compute_loss(current_response);
    train_scores(i_w,3) = alpha;
    
    
    MR = sum(((current_response>0)~=(labels>0)) .* wgt_samp)/...
        (length(labels) * mean(wgt_samp));
    
    train_scores_w(i_w,1) = 100*MR;

    
    if(i_w > 1)
       
        MR = sum((current_response_ctx1>0)~=(labels>0))/length(labels);
        
        fprintf('  Ctx1  Misclassif rate: %.2f | Loss: %f\n',100*MR,compute_loss(current_response_ctx1));

        train_scores_ctx1(i_w,1) = 100*MR;
      
        train_scores_ctx1(i_w,2) = compute_loss(current_response_ctx1);
        
        train_scores_ctx1(i_w,3) = alpha_ctx1;
        
        MR = sum(((current_response_ctx1>0)~=(labels>0)) .* wgt_samp)/...
            (length(labels) * mean(wgt_samp));
        
        train_scores_ctx1w(i_w,1) = 100*MR;
        
        save([params.codename '_train_scores_sav_w.mat'],'train_scores_w','train_scores_ctx1w');

        
        for i_img = 1 : length(data.train.gts.X)
            
            img_idx = find(samples_idx(:,1) == i_img);
            
            MR_img(i_w,i_img) = sum((current_response(img_idx)>0)~=(labels(img_idx)>0))/length(labels(img_idx));
            
            MR_img_ctx1(i_w,i_img) = sum((current_response_ctx1(img_idx)>0)~=(labels(img_idx)>0))/length(labels(img_idx));
            
        end
        
        save([params.codename '_MR_sav.mat'],'MR_img','MR_img_ctx1');
        
        save([params.codename '_train_scores_sav.mat'],'train_scores','train_scores_ctx1');
        
    end
    
    if(i_w > 0)
        
        
        % collect th colour value of the sample data
        
        wltree = weak_learners(i_w).reg_tree;
        
        for i_img = 1 : size(data.train.(ch).X,1)
            
            ch = params.ch_list{i_ch};
            
            X = data.train.(ch).X(i_img,data.train.(ch).idxs);
            
            gt_img = data.train.gts.X{i_img};
            
            clear I;
            
            for i_b = 1 : 3
                
                I(:,:,i_b) = X{i_b};
                
            end
            
            
            tStart = tic;
            
            
            % to avoid the computational burden, now apply the method only
            % on some sample points
                        
            samp_idx_img = samples_idx(samples_idx(:,1) == i_img,:);
            
            mask_img = data.train.masks.X{i_img};
            
            [mask_x,mask_y] = find(mask_img > 0);
            
            
            [score_image,leaf_image] = predict_img_wl_sample(X,params,weak_learners(i_w),[mask_x,mask_y]);
            
            t1 = toc(tStart);
            
            fprintf( ' Evaluating the  image %d took %d seconds \n', i_img, t1);
            
            
            % get to explore individual region recognised by individual leaf
            % node
            
            
            img_pg = size(leaf_image,1) * size(leaf_image,2);
            
            n_leaf = 0;
            
            ctx_c1 = [];
            
            n_samp_img = sum(leaf_image(:) > 0);
            
            for i_n = 1 : length(wltree)
                
                if(wltree(i_n).isLeaf)
                    
                    ctx_c1(end + 1) = i_n;
                    
                    n_leaf = n_leaf + 1;
                    
                    nv = wltree(i_n).value;
                    
                    n_idx = find(leaf_image == i_n);
                    
                    
                end
                
            end
            
            samp_idx_img = samples_idx(samples_idx(:,1) == i_img,2:3);
            
            ctx_ftrs1 = zeros(size(samp_idx_img,1),length(ctx_c1));
            
            samp_idx_img = sub2ind(size(leaf_image),samp_idx_img(:,1),samp_idx_img(:,2));
            
            for i_ctx = 1 : length(ctx_c1)
                
                tStart = tic;
                
                lf_ctx = (leaf_image == ctx_c1(i_ctx));
                
                if(sum(lf_ctx(:)))
                    
                    dist_ctx_map = bwdist(lf_ctx);
                    
                    ctx_ftrs1(:,i_ctx) = dist_ctx_map(samp_idx_img);
                    
                else
                    ctx_ftrs1(:,i_ctx) = 5000;
                    
                end
                
                t1 = toc(tStart);
                
                fprintf( ' Calculating the distance map of %d of the image %d took %d seconds \n', i_ctx, i_img, t1);
                
            end
            
            % set the pixels on the context as irrelevant to avoid the overfitting issues
            
            
            mask_img = data.train.masks.X{i_img};
            
            dist_ctx_map = bwdist(~mask_img);
            
            ctx_ftrs1(:,end + 1) = dist_ctx_map(samp_idx_img);
            
            known_mb = score_image > fig_thres;  
            
            known_mb = (bwdist(~mask_img) > 10) .* known_mb;
           
            if(sum(known_mb(:)))
              
                dist_ctx_map = bwdist(known_mb);
                
                ctx_ftrs1(:,end + 1) = dist_ctx_map(samp_idx_img);
            
            else
                
                ctx_ftrs1(:,end + 1) = 5000;
                
            end
            
            ctx_ftrs1_w{i_w,i_img} = ctx_ftrs1;
            
            
        end
        
    end
    
    wl_time = toc(t_wl);
    
    fprintf('  Learning WL %d took %f seconds\n------------------------------------------------\n\n',i_w,wl_time);

end

clf;
figure(1);
plot(1:params.wl_no,train_scores(:,1),'b')
legend('MR');
saveas(gcf,fullfile(params.results_dir,'MR_train_scores.jpg'),'jpg');
figure(2);
plot(1:params.wl_no,train_scores(:,2),'g');
legend('loss');
saveas(gcf,fullfile(params.results_dir,'LOSS_train_scores.jpg'),'jpg');
figure(3);
plot(1:params.wl_no,train_scores(:,3),'r');
legend('alpha');
saveas(gcf,fullfile(params.results_dir,'ALPHA_train_scores.jpg'),'jpg');

end
