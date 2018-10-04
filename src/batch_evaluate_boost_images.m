function score_images = batch_evaluate_boost_images(imgs,params,weak_learners,imgs_mask)

if(nargin < 4)
    
    for i_img = 1 : size(imgs,1)
        
        imgs_mask{i_img} = ones(size(imgs{1}(:,:,1)));
        
    end
    
end

border_img = zeros(size(imgs{1}(:,:,1)));

border_img(1:params.border_skip_size,:) = 1;

border_img(:,1:params.border_skip_size) = 1;

border_img((end-params.border_skip_size):end,:) = 1;

border_img(:,(end-params.border_skip_size):end) = 1;

border_img = ~border_img;

roi_idxs = [];


for i_img = 1 : size(imgs,1)
    
    img_mask = imgs_mask{i_img} .* border_img;
    
    mask_idx = find(img_mask > 0);
    
    roi_idx = ones(length(mask_idx),3) * i_img;
    
    [roi_idx(:,2),roi_idx(:,3)] = ind2sub(size(imgs{1}(:,:,1)),mask_idx);
    
    roi_idxs = cat(1,roi_idxs,roi_idx);
    
end

wl_no = length(weak_learners);

score = zeros(size(roi_idxs,1),1);

sub_ch_no = size(imgs,2);


for i_w = 1:wl_no
    fprintf('    WL %d/%d\n',i_w,wl_no);
    
    wl = weak_learners(i_w);
    
    features = zeros(size(roi_idxs,1),length(wl.kernels));
    
    for i_ch = 1:params.ch_no
        ch = params.ch_list{i_ch};
        
        % Prepare the data
        wl = weak_learners(i_w);
       
        for i_s = 1:sub_ch_no
            
            idxs = find(cellfun(@(x)(x.ch_no==i_ch && x.sub_ch_no==i_s),...
                weak_learners(i_w).kernel_params));
            
            if (~isempty(idxs))
            
                features(:,idxs) = mexEvaluateKernels(imgs(:,i_s),roi_idxs(:,1:3),...
                    params.sample_size,weak_learners(i_w).kernels(idxs),weak_learners(i_w).kernel_params(idxs));
            
            end
            
            
        end

    end
    
    new_score = LDARegStumpPredict(wl.reg_tree,single(features));
    
    score = score + wl.alpha * new_score;
    
end

s_m = min(score(:));

s_M = max(score(:));

score = (score - s_m)/(s_M - s_m);


for i_img = 1 : size(imgs,1)
    
    score_image =  zeros(size(imgs{i_img}(:,:,1)));
    
    score_img = score(roi_idxs(:,1) == i_img);
    
    roi_idx_img = roi_idxs(roi_idxs(:,1) == i_img,2:3);
    
    roi_idx_img = sub2ind(size(score_image),roi_idx_img(:,1),roi_idx_img(:,2));
    
    score_image(roi_idx_img) = score_img;
    
    score_images{i_img} = score_image;
    
end

