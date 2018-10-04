%
% samples_idx(:,1) => sample image no
% samples_idx(:,2) => sample row
% samples_idx(:,3) => sample column
% samples_idx(:,4) => sample label (-1/+1)
function [weak_learners,train_scores] = train_boost_3D_CLRG_combined(params,data,ftrs,wgt,samples_idx,LTClassifier)

% Train a KernelBoost classifier on the given samples

% the classifier combine the added discriptor


% allows an additional weight, i.e. assign the weight according to the
% number of pixels consisted in the seed

% works on 3D iamge, in accordance with the paper, combine with a CLRG tree
% classifier


%
%  authors: Carlos Becker, Roberto Rigamonti, CVLab EPFL
%  e-mail: name <dot> surname <at> epfl <dot> ch
%  web: http://cvlab.epfl.ch/
%  date: February 2014

samples_no = size(samples_idx,1);

weak_learners(params.wl_no).alpha = 0;
labels = samples_idx(:,5);
samples_idx = samples_idx(:,1:4);

samples_idx(:,2:3) = samples_idx(:,2:3) + params.border_size;


current_response = zeros(samples_no,1);

[compute_wi,compute_ri,compute_loss,compute_indiv_loss,compute_2nd_deriv,mex_loss_type] = select_fncts(params,labels);

W = compute_wi(current_response);

W = W .* wgt;


R = compute_ri(current_response);

train_scores = zeros(params.wl_no,4);

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
        
        X = expand_img(data.train.(ch).X,params);
        
        
%         for i_img = 1 : length(data.train.(ch).X)
%             
%             X{i_img,1} = sum(data.train.(ch).X{i_img},3);
%             
%         end
        
        
        %X = data.train.(ch).X(:,data.train.(ch).idxs);
        X_idxs = 1;
        
        sub_ch_no = data.train.(ch).sub_ch_no;
        
        sub_ch_no = 1;
        
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
    
    % add the histogram discriptor
    
    ftrs1 = ftrs(T2_idx,:);

    features = [features,ftrs1];
    
    
    fprintf('    Training regression tree on learned features...\n');
    t_tr = tic;
    reg_tree = LDARegStumpTrain(single(features),R(T2_idx),W(T2_idx)/sum(W(T2_idx)),uint32(params.tree_depth));
    time_tr = toc(t_tr);
    fprintf('    Done! (took %f seconds)\n',time_tr);
    
    fprintf('    Removing useless kernels...\n');
    [weak_learners(i_w).kernels,weak_learners(i_w).kernel_params,weak_learners(i_w).reg_tree...
        ] = remove_useless_filters_ftrs(reg_tree,kernels,kernel_params);
    
        
    t_ev = tic;
    fprintf('    Evaluating the learned kernels on the whole training set...\n');
    features = zeros(length(labels),length(weak_learners(i_w).kernels));
    for i_ch = 1:params.ch_no
        ch = params.ch_list{i_ch};
        sub_ch_no = data.train.(ch).sub_ch_no;
        
        sub_ch_no = 1;
        
       % X = data.train.(ch).X(:,data.train.(ch).idxs);
        
        X = expand_img(data.train.(ch).X,params);
        
        for i_s = 1:sub_ch_no
            idxs = find(cellfun(@(x)(x.ch_no==i_ch && x.sub_ch_no==i_s),weak_learners(i_w).kernel_params));
%             idxs = 1;
            
            if (~isempty(idxs))
                features(:,idxs) = mexEvaluateKernels(X(:,i_s),samples_idx(:,1:3),params.sample_size,weak_learners(i_w).kernels(idxs),weak_learners(i_w).kernel_params(idxs));
            end
        end
    end
    ev_time = toc(t_ev);
    fprintf('      Evaluation completed in %f seconds\n',ev_time);
    

    % add the hd feature 
    
    
    features = [features,ftrs];    
        
    fprintf('    Performing prediction on the whole training set...\n');
    t_pr = tic;
    cached_responses = LDARegStumpPredict(weak_learners(i_w).reg_tree,single(features));
    time_pr = toc(t_pr);
    fprintf('    Prediction finished, took %f seconds\n',time_pr);
    clear features;
    
    fprintf('    Finding alpha through line search...\n');
    t_alp = tic;
    alpha = mexLineSearch(current_response,cached_responses,labels,mex_loss_type);
    time_alp = toc(t_alp);
    fprintf('    Good alpha found (alpha=%f), took %f seconds\n',alpha,time_alp);
    alpha = alpha * params.shrinkage_factor;
    
    current_response = current_response + alpha*cached_responses;
    
    W = compute_wi(current_response);
    
    % incorperate the weight
    W = W .* wgt;
    
    R = compute_ri(current_response);
    
    weak_learners(i_w).alpha = alpha;
    
    MR = sum((current_response>0)~=(labels>0))/length(labels);
    fprintf('    Misclassif rate: %.2f | Loss: %f\n',100*MR,compute_loss(current_response));
    
    WMR = sum(((current_response>0)~=(labels>0)) .* wgt )/sum(wgt);
    fprintf('  Weighted  Misclassif rate: %.2f | Loss: %f\n',100*WMR,compute_loss(current_response));    
    
    train_scores(i_w,1) = 100*MR;
    train_scores(i_w,2) = compute_loss(current_response);
    train_scores(i_w,3) = alpha;
    train_scores(i_w,4) = 100 * WMR;
    
    
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

figure(4);
plot(1:params.wl_no,train_scores(:,4),'r');
legend('WMR');
saveas(gcf,fullfile(params.results_dir,'Weighted_MR_train_scores.jpg'),'jpg');

end
