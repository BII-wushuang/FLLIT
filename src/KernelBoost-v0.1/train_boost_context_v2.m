%
% samples_idx(:,1) => sample image no
% samples_idx(:,2) => sample row
% samples_idx(:,3) => sample column
% samples_idx(:,4) => sample label (-1/+1)
function [weak_learners] = train_boost_context_v2(params,data,samples_idx)

% Train a KernelBoost classifier on the given samples
%
%  authors: Carlos Becker, Roberto Rigamonti, CVLab EPFL
%  e-mail: name <dot> surname <at> epfl <dot> ch
%  web: http://cvlab.epfl.ch/
%  date: February 2014

samples_no = size(samples_idx,1);

% weak_learners(params.wl_no).alpha = 0;
% 
% weak_learners_ctx1(params.wl_no).alpha = 0;
% 
% weak_learners_ctx2(params.wl_no).alpha = 0;




labels = samples_idx(:,4);
samples_idx = samples_idx(:,1:3);
current_response = zeros(samples_no,1);

[compute_wi,compute_ri,compute_loss,compute_indiv_loss,compute_2nd_deriv,mex_loss_type] = select_fncts(params,labels);

W = compute_wi(current_response);
R = compute_ri(current_response);

current_response_ctx1 = current_response;

current_response_ctx2 = current_response;

W_ctx1 = W;

R_ctx1 = R;


W_ctx2 = W;

R_ctx2 = R;


train_scores = zeros(params.wl_no,3);

train_scores_ctx1 = zeros(params.wl_no,3);

train_scores_ctx2 = zeros(params.wl_no,3);

n_km = 30;


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
        

        
        %save('tmp_ftrs.mat');
        
        nf = size(features,1);
        
        ncf_g = size(ctx_ftrs1,2);
        
        ctx_glb_all = zeros(size(samples_idx,1),ncf_g);
        
        ctx_loc_all = zeros(size(samples_idx,1),n_ftrs2);
        
        save('tmp_ftrs.mat');
        
        
        for i_img = 1 : size(data.train.(ch).X,1)
            
            ctx_glb_all(samples_idx(:,1) == i_img,:) = ctx_ftrs1_w{i_w - 1,i_img};
            
            if(isempty(ctx_ftrs2_w{i_w - 1,i_img}))
            
                save('empty_flag_sav.mat', ctx_ftrs2_w);
                
                %ctx_loc_all(samples_idx(:,1) == i_img,:) = ones(sum(samples_idx(:,1) == i_img),n_km) * 5000;
                
            else
                
                ctx_loc_all(samples_idx(:,1) == i_img,:) = ctx_ftrs2_w{i_w - 1,i_img};
            
            end
            
            ctx_ftrs1_w{i_w - 1,i_img} = [];
            
            ctx_ftrs2_w{i_w - 1,i_img} = [];
            
        end
        
        ctx_glb = ctx_glb_all(T2_idx,:);
        
        ctx_loc = ctx_loc_all(T2_idx,:);
        
        
        fprintf('    Training regression tree on learned features combined with context 1...\n');
        t_tr = tic;
        reg_tree_g = LDARegStumpTrain(single([features,...
            ctx_glb]),R_ctx1(T2_idx),W_ctx1(T2_idx)...
            /sum(W_ctx1(T2_idx)),uint32(params.tree_depth));
        time_tr = toc(t_tr);
        fprintf('    Done! (took %f seconds)\n',time_tr);
        
        
        %     t_tr = tic;
        %     reg_tree_l = LDARegStumpTrain(single([features,ctx_loc]),R(T2_idx),W(T2_idx)/sum(W(T2_idx)),uint32(params.tree_depth));
        %     time_tr = toc(t_tr);
        %     fprintf('    Done! (took %f seconds)\n',time_tr);
        
        
        
        %     for i_img = 1 : size(data.train.(ch).X,1)
        %
        %         T2_img = samples_idx(T2_idx,1) == i_img;
        %
        %         ctx_loc_img = ctx_loc(T2_img,:);
        %
        %         ctx_glb_img = ctx_glb(T2_img,:);
        %
        %         fprintf('    Training regression tree %d on learned features combined with context 1 2...\n',i_img);
        %
        %         t_tr = tic;
        %         reg_tree_l{i_img} = LDARegStumpTrain(single([features(T2_img,:),...
        %             ctx_glb_img,ctx_loc_img]),R(T2_idx(T2_img)),W(T2_idx(T2_img)) ...
        %             /sum(W(T2_idx((T2_img)))),...
        %             uint32(params.tree_depth));
        %         time_tr = toc(t_tr);
        %         fprintf('    Done! (took %f seconds)\n',time_tr);
        %
        %     end
        
        fprintf('    Training regression tree on learned features combined with context 1 2...\n');
        
        t_tr = tic;
        
        reg_tree_l = LDARegStumpTrain(single([features,...
            ctx_glb,ctx_loc]),R_ctx2(T2_idx),W_ctx2(T2_idx) ...
            /sum(W_ctx2(T2_idx)),uint32(params.tree_depth));
        
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
        
        [weak_learners_ctx2(i_w).kernels,weak_learners_ctx2(i_w).kernel_params,weak_learners_ctx2(i_w).reg_tree,...
            weak_learners_ctx2(i_w).ctx_list]...
            = remove_useless_filters_ctx(reg_tree_l,kernels,kernel_params);
        
        ctx2_list = weak_learners_ctx2(i_w).ctx_list;
        
        ctx2_list(ctx2_list < ncf_g + 1) = [];
        
        ctx2_list = ctx2_list - ncf_g;
        
        tmp_kmc = [];
        
        tmp_ltree = [];
        
        for i_km = 1 : length(ctx2_list)
           
            tmp_kmc(i_km,:) = kmc_ctx_c2(ctx2_list(i_km),:);
            
            tmp_ltree(i_km) = tleaf_ctx_c2(ctx2_list(i_km));
            
        end
        
        weak_learners_ctx2(i_w).kmc = tmp_kmc;
        
        weak_learners_ctx2(i_w).ltree = tmp_ltree;
        
       % ctx2_list = ctx2_list - ncf_g;
        
        weak_learners_ctx2(i_w).ctx2_list = ctx2_list;
        
        save('weak_learners_sav.mat','weak_learners_ctx2','weak_learners_ctx1','weak_learners');
            
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
        
        ctx_glb_all = [ctx_glb_all, ctx_loc_all];
        
        clear ctx_loc_all;
        
        features_ctx2 = ctx_glb_all(:,weak_learners_ctx2(i_w).ctx_list);
        
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
        
        
        t_ev = tic;
        fprintf('    Evaluating the learned kernels on the whole training set...\n');
        features = zeros(length(labels),length(weak_learners_ctx2(i_w).kernels));
        for i_ch = 1:params.ch_no
            ch = params.ch_list{i_ch};
            sub_ch_no = data.train.(ch).sub_ch_no;
            
            X = data.train.(ch).X(:,data.train.(ch).idxs);
            for i_s = 1:sub_ch_no
                idxs = find(cellfun(@(x)(x.ch_no==i_ch && x.sub_ch_no==i_s),...
                    weak_learners_ctx2(i_w).kernel_params));
                if (~isempty(idxs))
                    features(:,idxs) = mexEvaluateKernels(X(:,i_s),...
                        samples_idx(:,1:3),params.sample_size,...
                        weak_learners_ctx2(i_w).kernels(idxs),...
                        weak_learners_ctx2(i_w).kernel_params(idxs));
                end
            end
        end
        ev_time = toc(t_ev);
        fprintf('      Evaluation completed in %f seconds\n',ev_time);        
        
        

        fprintf('    Performing ctx2 prediction on the whole training set...\n');
        
        t_pr = tic;
        cached_responses_ctx2 = LDARegStumpPredict(...
            weak_learners_ctx2(i_w).reg_tree,single([features, features_ctx2]));
        time_pr = toc(t_pr);
        fprintf('    Prediction finished, took %f seconds\n',time_pr);        
        
        clear features;
        
        clear features_ctx1 features_ctx2 ctx_glb_all;
        
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
        
        
        fprintf('    Finding alpha 2 through line search...\n');
        t_alp = tic;
        alpha_ctx2 = mexLineSearch(current_response_ctx2,cached_responses_ctx2,labels,mex_loss_type);
        time_alp = toc(t_alp);
        fprintf('    Good alpha found (alpha=%f), took %f seconds\n',alpha_ctx2,time_alp);
        
        
        
    end
    
    
    alpha = alpha * params.shrinkage_factor;
    
    current_response = current_response + alpha*cached_responses;
    
    W = compute_wi(current_response);
    R = compute_ri(current_response);
    
    weak_learners(i_w).alpha = alpha;
    
    
    if(i_w > 1)
        
        alpha_ctx1 = alpha_ctx1 * params.shrinkage_factor;
        
        current_response_ctx1 = current_response_ctx1 + alpha_ctx1*cached_responses_ctx1;
        
        W_ctx1 = compute_wi(current_response_ctx1);
        
        R_ctx1 = compute_ri(current_response_ctx1);
        
        weak_learners_ctx1(i_w).alpha = alpha_ctx1;
       
        
        alpha_ctx2 = alpha_ctx2 * params.shrinkage_factor;
        
        current_response_ctx2 = current_response_ctx2 + alpha_ctx2*cached_responses_ctx2;
        
        W_ctx2 = compute_wi(current_response_ctx2);
        
        R_ctx2 = compute_ri(current_response_ctx2);
        
        weak_learners_ctx2(i_w).alpha = alpha_ctx2;       
        
        
        save('weak_learners_sav.mat','weak_learners','weak_learners_ctx1','weak_learners_ctx1')
        
        
    end
    
    MR = sum((current_response>0)~=(labels>0))/length(labels);
    fprintf('    Misclassif rate: %.2f | Loss: %f\n',100*MR,compute_loss(current_response));
    
    train_scores(i_w,1) = 100*MR;
    train_scores(i_w,2) = compute_loss(current_response);
    train_scores(i_w,3) = alpha;
    
    
    if(i_w > 1)
       
        MR = sum((current_response_ctx1>0)~=(labels>0))/length(labels);
        
        fprintf('  Ctx1  Misclassif rate: %.2f | Loss: %f\n',100*MR,compute_loss(current_response_ctx1));

        train_scores_ctx1(i_w,1) = 100*MR;
      
        train_scores_ctx1(i_w,2) = compute_loss(current_response_ctx1);
        
        train_scores_ctx1(i_w,3) = alpha_ctx1;
        
        
        
        MR = sum((current_response_ctx2>0)~=(labels>0))/length(labels);
        
        fprintf('  Ctx2  Misclassif rate: %.2f | Loss: %f\n',100*MR,compute_loss(current_response_ctx2));
      
        train_scores_ctx2(i_w,1) = 100*MR;
      
        train_scores_ctx2(i_w,2) = compute_loss(current_response_ctx2);
        
        train_scores_ctx2(i_w,3) = alpha_ctx2;        
        
        
    end
       
    save('train_scores_sav.mat','train_scores','train_scores_ctx1','train_scores_ctx2');
    
    if(i_w > 0)
        
        
        % collect th colour value of the sample data
        
        samp_rgb = zeros(size(samples_idx,1),3);
        
        for i_img = 1 : size(data.train.(ch).X,1)
            
            X = data.train.(ch).X(i_img,data.train.(ch).idxs);
            
            clear I;
            
            for i_b = 1 : 3
                
                I(:,:,i_b) = X{i_b};
                
            end
            
            
            samp_idx_img = samples_idx(samples_idx(:,1) == i_img,:);
            
            posi_samp = sub2ind(size(I(:,:,1)),samp_idx_img(:,2),samp_idx_img(:,3));
            
            I = reshape(I,[],3);
            
            samp_rgb(samples_idx == i_img,:) = I(posi_samp,:);
            
        end
        
        
        wltree = weak_learners(i_w).reg_tree;
        
        [~,samp_lf] = predict_idx_wl(data,params,samples_idx,weak_learners(i_w));
        
        lf_hist = histc(samp_lf,1:length(wltree));
        
        lf_hist = lf_hist(1:length(wltree));
        
        lf_hist = lf_hist / sum(lf_hist);
        
        [lfn,lf_id] = sort(lf_hist,'descend');
        
        lf_id(lfn < 0.3) = [];
        
        lfn(lfn < 0.3) = [];
        

        
        %     for i_lf = 1 : length(lf_id)
        %
        %         [~,ctx_km_c{i_lf}] = kmeans(samp_rgb(samp_lf == lf_id(i_lf),:),n_km,'EmptyAction','singleton');
        %
        %     end
        %
        
        for i_n = 1 : length(wltree)
            
            bkg_ftrs2{i_n} = [];
            
            
            
            %bkg_idxs{i_n} = [];
            
        end
        
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
            
            [score_image,leaf_image] = predict_img_wl_sample(X,params,weak_learners(i_w),samp_idx_img(:,2:3));
            
            t1 = toc(tStart);
            
            fprintf( ' Evaluating the  image %d took %d seconds \n', i_img, t1);
            
            
            % get to explore individual region recognised by individual leaf
            % node
            
            
            img_pg = size(leaf_image,1) * size(leaf_image,2);
            
            %context_img = zeros(size(leaf_image));
            
            %   cxt_idx = 1;
            
            
            n_leaf = 0;
            
            ctx_c1 = [];
            
            ctx_c2 = [];
            
            
            n_samp_img = sum(leaf_image(:) > 0);
            
            for i_n = 1 : length(wltree)
                
                if(wltree(i_n).isLeaf)
                    
                    ctx_c1(end + 1) = i_n;
                    
                    n_leaf = n_leaf + 1;
                    
                    nv = wltree(i_n).value;
                    
                    n_idx = find(leaf_image == i_n);
                    
                    
                    if(length(n_idx) / n_samp_img > 0.2)
                        
                        I_2D = reshape(I,[],3);
                        
                        n_idx = downsample(n_idx,200);
                        
                        bkg_ftrs2{i_n} = [ bkg_ftrs2{i_n}; I_2D(n_idx,:)];
                        
                        
                    end
                    
                end
                
            end
            
            % contx_list{i_w,i_img} = ctx_c2;
            
            samp_idx_img = samples_idx(samples_idx(:,1) == i_img,2:3);
            
            ctx_ftrs1 = zeros(size(samp_idx_img,1),length(ctx_c1));
            
            % ctx_ftrs2 = zeros(size(samp_idx_img,1),length(ctx_c2));
            
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
            
            ctx_ftrs1(ctx_ftrs1 < 0.1) = 5000;
            
            ctx_ftrs1_w{i_w,i_img} = ctx_ftrs1;
            
        end
        
        
        n_km = 30;
        
        ctx_c2 = [];
        
        clear ftrs_kmc_c2;
        
        kmc_ctx_c2 = [];
        
        tleaf_ctx_c2 = [];
        
        for i_n = 1 : length(wltree)
            
            if(~isempty(bkg_ftrs2{i_n}))
                
                tStart = tic;
                
                [~,ftrs_kmc_tmp] = kmeans(bkg_ftrs2{i_n},n_km,'EmptyAction','singleton');
                
                t1 = toc(tStart);
                
                fprintf('Clustering the context ftrs 2 took %d seconds', t1);
                
                
                ctx_c2(end + 1 : end + n_km) = (100 * i_n) + (1:n_km);
                
                ftrs_kmc{i_w,i_n} = ftrs_kmc_tmp;
                
                kmc_ctx_c2(end + 1: end + n_km,:) = ftrs_kmc_tmp;
                
                tleaf_ctx_c2(end + 1: end + n_km,:) = i_n;
                
            else
                
                ftrs_kmc{i_w,i_n} = {};
                
            end
            
        end
        
        n_ftrs2 = length(ctx_c2);
        
        
        for i_img = 1 : size(data.train.(ch).X,1)
            
            
            ctx_ftrs2_w{i_w,i_img} = [];
            
            tStart = tic;
            
            
            % to avoid the computational burden, now apply the method only
            % on some sample points
                        
            samp_idx_img = samples_idx(samples_idx(:,1) == i_img,:);
            
            [score_image,leaf_image] = predict_img_wl_sample(X,params,weak_learners(i_w),samp_idx_img(:,2:3));
            
            t1 = toc(tStart);
            
            fprintf( ' Evaluating the  image %d took %d seconds \n', i_img, t1);
            
            ftrs2_img = zeros(size(leaf_image));
            
            samp_idx_img = samples_idx(samples_idx(:,1) == i_img,2:3);
            
            samp_idx_img = sub2ind(size(leaf_image),samp_idx_img(:,1),samp_idx_img(:,2));
            
            ctx_ftrs2 = zeros(size(samp_idx_img,1),length(ctx_c2));
            
            for i_n = 1 : length(wltree)
                
                if(~isempty(ftrs_kmc{i_w,i_n}))
                    
                    ftrs2_img_tmp = ftrs_c2_img(X,ftrs_kmc{i_w,i_n});
                    
                    ftrs2_img = (i_n * 100 + ftrs2_img_tmp) .* (leaf_image == i_n);
                    
                    
                end
                
            end
            
            
            for i_ctx = 1 : length(ctx_c2)
                
                tStart = tic;
            
                lf_ctx = (ftrs2_img == ctx_c2(i_ctx));
                
                if(sum(lf_ctx(:)))
                    
                    dist_ctx_map = bwdist(lf_ctx);
                    
                    ctx_ftrs2(:,i_ctx) = dist_ctx_map(samp_idx_img);
                    
                else
                    
                    ctx_ftrs2(:,i_ctx) = 5000;
                    
                end
                
                t1 = toc(tStart);
                
                fprintf( ' Calculating the distance map of %d of the image %d took %d seconds \n', i_ctx, i_img, t1);
                
            end
            
            ctx_ftrs2(ctx_ftrs2 < 0.1) = 5000;
            
            ctx_ftrs2_w{i_w,i_img} = ctx_ftrs2;
            
            %  clear ctx_ftrs2;
            
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
