function [] = test_classifiers_output_imgs(params,data,weak_learners)

% Test KernelBoost classifier on data
%
%  authors: Carlos Becker, Roberto Rigamonti, CVLab EPFL
%  e-mail: name <dot> surname <at> epfl <dot> ch
%  web: http://cvlab.epfl.ch/
%  date: February 2014

parfor i_img = 1:params.gts.test.imgs_no
    t_cl = tic;
    fprintf('--\n  Predicting img %d/%d\n',i_img,params.gts.test.imgs_no); %#ok<*PFBNS>
    tmp_image = data.test.gts.X{i_img,1};
    
    score_image = zeros(size(tmp_image));
    wl_no = length(weak_learners);
    
    [~,xx,~] = fileparts(params.gts.test.list{i_img});
    score_filename = fullfile(params.test_subdir_path,sprintf('%s.nrrd',xx));
    
    if (~exist(score_filename,'file'))
        for i_w = 1:wl_no
            fprintf('    WL %d/%d\n',i_w,wl_no);
            wl = weak_learners(i_w);
            features = zeros(numel(tmp_image),length(wl.kernels));
            for i_ch = 1:params.ch_no
                ch = params.ch_list{i_ch};
                
                % Prepare the data
                wl = weak_learners(i_w);
                idxs = find(cellfun(@(x)(strcmp(x.ch_name,ch)),wl.kernel_params));
                if (length(idxs)>0)
                    tmp_kp = wl.kernel_params{idxs(1)};
                    X = data.test.(ch).X(i_img,tmp_kp.data_idxs);
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
            score_image = score_image+reshape(wl.alpha*new_score,[size(tmp_image,1),size(tmp_image,2)]);
        end
        score_image = 1./(1+exp(-2.*score_image));
    else
        score_image = nrrdLoad(score_filename)';
    end
    
    reduced_mask = data.test.masks.X{i_img}(params.border_size+1:end-params.border_size,params.border_size+1:end-params.border_size);
    score_image = score_image(params.border_size+1:end-params.border_size,params.border_size+1:end-params.border_size);
    score_image(reduced_mask==0) = min(score_image(reduced_mask==1));
    s_m = min(score_image(:));
    s_M = max(score_image(:));
    score_image = (score_image-s_m)/(s_M-s_m);
    nrrdSave(score_filename,permute(score_image,[2,1]));
    
    imwrite();
    
    time_tcl = toc(t_cl);
    fprintf('  Done, took %f seconds\n',time_tcl);
end

end
