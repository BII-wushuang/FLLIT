function new_score = KB_evaluate(X,wl,params)

features = zeros(numel(X{1}(:)),length(wl.kernels));

for i_ch = 1:params.ch_no
    ch = params.ch_list{i_ch};
    
    % Prepare the data
    idxs = find(cellfun(@(x)(strcmp(x.ch_name,ch)),wl.kernel_params));
    if (length(idxs)>0)
        tmp_kp = wl.kernel_params{idxs(1)};
        X = X(tmp_kp.data_idxs);
    else
        X = [];
    end
    ch_wl = wl;
    ch_wl.kernels = ch_wl.kernels(idxs);
    ch_wl.kernel_params = ch_wl.kernel_params(idxs);
    
    for i_k = 1:length(ch_wl.kernels)
        kernel_params = ch_wl.kernel_params{i_k};
        k_s = kernel_params.filter_size;
        r_k = kernel_params.start_row;
        c_k = kernel_params.start_col;
        kernel = zeros(params.sample_size(1),params.sample_size(2));
        kernel(r_k:r_k+k_s-1,c_k:c_k+k_s-1) = reshape(ch_wl.kernels{i_k},[k_s,k_s]);
        response = imfilter(X{ch_wl.kernel_params{i_k}.sub_ch_no},kernel,'same');
        features(:,idxs(i_k)) = response(:);
    end
end

new_score = LDARegStumpPredict(wl.reg_tree,single(features));
           % score_image = score_image+reshape(wl.alpha*new_score,[size(tmp_image,1),size(tmp_image,2)]);