%Classifier for leg segmentation
function Segmentation (data_dir,score_thres,foreground_thres,load_wl)
%% locate the image folders and the output folders
if (nargin < 1)
    load_wl = 1;
    score_thres = 0.65;
    foreground_thres = 0.1;
    data_dir = uigetdir('./Data');
    addpath(genpath('./KernelBoost-v0.1/'));
end

clear weak_learners;
%addpath(genpath('./KernelBoost-v0.1/'));

%fprintf('Processing the folder: %s \n',data_dir);

pos_bs = strfind(data_dir,'Data');
sub_dir = data_dir(pos_bs(end)+length('Data'):length(data_dir));
data_dir = [pwd '/Data' sub_dir '/'];
output_dir = [pwd '/Results/SegmentedImages' sub_dir '/'];
if(~exist(output_dir))
    mkdir(output_dir);
end

img_list = load_img_list(data_dir);
I = imread([data_dir img_list(1).name]);

%% Call the background 
if(~exist([data_dir 'Background/Background.png']))
    [~,~,ref_img,~] = video2background(data_dir, sub_dir);
    if(~exist([data_dir 'Background']))
        mkdir([data_dir 'Background'])
    end
    imwrite(uint8(ref_img),[data_dir 'Background/Background.png'], 'png');
else
    ref_img = imread([data_dir 'Background/Background.png']);
end

imshow(uint8(ref_img));
pause(2);

for i = 1 : length(img_list)
    I = imread([data_dir img_list(i).name]);
    I = double(I);

    roi_img = (max(ref_img - I(:,:,1),0) ./ ref_img) > foreground_thres;
    roi_img = bwareaopen(roi_img, 200);
    
    img_output = repmat(I/255,[1 1 3]);
    img_output(:,:,3) = img_output(:,:,3) + roi_img;
    imshow(img_output);
    pause(0.01);
    imwrite(roi_img,[output_dir 'roi_' num2str(i) '.png'],'png');
end

%% Another section

ref_img = imread([data_dir 'Background/Background.png']);
ref_img = double(ref_img);
ref_img = padarray(ref_img,[20 20],'replicate');

params = setup_config_L('Drive');
%params = setup_lists(params);
%params = setup_directories_L(params);
params.border_skip_size = 20;
params.pos_samples_no = 30000;
params.neg_samples_no = 30000;
params.sample_size = [41 41]';

if (load_wl)
    fprintf('Please select the trained classifier to be used.\n');
    [wl_file, wl_dir] = uigetfile([pwd '/Results/Classifiers/','*.mat']);
    load([wl_dir wl_file]);
else    
    if(~exist([pwd '/Results/Classifiers' sub_dir '_Classifier.mat']))
        sample_ratio = 20;
        tot_idx.pos = 0;
        tot_idx.neg = 0; 
        data = [];

        % Collect the positive and negative samples
        for i_img = 1 : floor(length(img_list) / sample_ratio)

            I = imread([data_dir img_list(i_img * sample_ratio).name]); 
            I = double(I);
            I = padarray(I,[20 20],'replicate');

            [pos_img,neg_img_body,neg_img_bkg] = leg_segment(I,ref_img,foreground_thres);
            show_img_output = repmat(I/255,[1 1 3]);
            show_img_output(:,:,1) = show_img_output(:,:,1) + pos_img;
            show_img_output(:,:,3) = show_img_output(:,:,3) + neg_img_body;
            imshow(show_img_output,[]);

            data.train.imgs.X{i_img,1} = I / 255;
            data.train.imgs.X{i_img,2} = (ref_img - I) / 255;

            gt_img = zeros(size(I));  
            gt_img(pos_img > 0) = 1;    
            gt_img(neg_img_body > 0) = -1;      
            gt_img(neg_img_bkg > 0) = -1;
            data.train.gts.X{i_img} = gt_img;

            % ------------------------------

            border_img = zeros(size(I));
            border_img(1:params.border_skip_size,:) = 1;
            border_img(:,1:params.border_skip_size) = 1;
            border_img((end-params.border_skip_size):end,:) = 1;
            border_img(:,(end-params.border_skip_size):end) = 1;

            pos_sampling = find((pos_img) & (~border_img) & (~neg_img_body)); 
            neg_body_sampling = find(neg_img_body & (~border_img) & (~pos_img));
            neg_bkg_sampling = find(neg_img_bkg & (~border_img));

            % sample counts
            imgs_no = floor(length(img_list) / sample_ratio);
            npos_img = ceil(params.pos_samples_no / imgs_no);
            nneg_img = ceil(params.neg_samples_no / imgs_no);

            % getting a random subset?
            neg_sampling = neg_body_sampling(randi(length(neg_body_sampling),[nneg_img,1]));
            neg_sampling = [neg_sampling; neg_bkg_sampling(randi(length(neg_bkg_sampling),[nneg_img,1]))];

            npos = length(pos_sampling);        
            nneg = length(neg_sampling);

            good_sampling_idx{i_img}.pos = pos_sampling(randi(npos,[max(npos_img,nneg_img),1]));        
            good_sampling_idx{i_img}.neg = neg_sampling(randi(nneg,[max(npos_img,nneg_img),1]));
            good_sampling_idx{i_img}.pos_val = pos_img(good_sampling_idx{i_img}.pos);        
            good_sampling_idx{i_img}.neg_val = pos_img(good_sampling_idx{i_img}.neg);

            tot_idx.pos = tot_idx.pos+length(good_sampling_idx{i_img}.pos);        
            tot_idx.neg = tot_idx.neg+length(good_sampling_idx{i_img}.neg);
        end
    
        data.train.imgs.idxs = 1:2;
        data.train.imgs.sub_ch_no = 2;
        samples_idx = [];

        for i_img = 1 : floor(length(img_list) / sample_ratio)

            I = imread([data_dir img_list(i_img * sample_ratio).name]);
            I = padarray(I,[20 20],'replicate');

            samp_idx = [good_sampling_idx{i_img}.pos ; good_sampling_idx{i_img}.neg];

            samp_idx_2D = zeros(length(samp_idx),2);
            [samp_idx_2D(:,1),samp_idx_2D(:,2)] = ind2sub(size(I),samp_idx);

            labels = [good_sampling_idx{i_img}.pos_val ; good_sampling_idx{i_img}.neg_val];
            labels = double(labels);
            labels(labels < 0.5) = -1;

            samples_idx_img = zeros(length(samp_idx),4);
            samples_idx_img(:,1) = i_img;
            samples_idx_img(:,2:3) = samp_idx_2D;
            samples_idx_img(:,4) = labels;

            samples_idx = cat(1,samples_idx,samples_idx_img);

        end

        samples_idx = samples_idx(1 : (params.pos_samples_no + params.neg_samples_no),:);

        params.pos_samples_no = sum(samples_idx(:,end)==1);   
        params.neg_samples_no = sum(samples_idx(:,end)==-1);  

        params.T1_size = round((params.pos_samples_no+params.neg_samples_no)/3);
        params.T2_size = round((params.pos_samples_no+params.neg_samples_no)/3);

        params.pos_to_sample_no = params.T1_size;
        params.neg_to_sample_no = params.T1_size;
        params.wl_no = 100;

        fprintf('Training the classifier over a random subset of the images.\n');
        % Train the classifier here
        weak_learners = train_boost_general(params,data,samples_idx);
        sub_dir_pos = strfind(sub_dir,'/');
        if(~exist(['./Results/Classifiers/' sub_dir(1:sub_dir_pos(end))]))
            mkdir(['./Results/Classifiers/' sub_dir(1:sub_dir_pos(end))]);
        end
        save ([pwd '/Results/Classifiers/' sub_dir '_Classifier.mat'],'weak_learners');
    else
        load([pwd '/Results/Classifiers/' sub_dir '_Classifier.mat']);
    end    
end

fprintf('Training completed, now applying the classifier to all images.\n');
fprintf('Output directory: %s\n', output_dir);

% Applying the classifier
sample_ratio = 1;
sec_no = 100;

for i_sec = 1 : ceil(length(img_list) / sample_ratio / sec_no)
    X = [];
    if (i_sec < ceil(length(img_list) / sample_ratio / sec_no))
        imgs_sec = 1 + (i_sec - 1) * sec_no : sample_ratio :(i_sec) * sec_no;
    else
        imgs_sec = 1 + (i_sec - 1) * sec_no : sample_ratio : length(img_list);
    end
    for i_img = 1 : length(imgs_sec)
        I = imread([data_dir img_list(imgs_sec(i_img)).name]);
        I = double(I);
        I = padarray(I,[20 20],'replicate');
        
        %neg_img{i_img} = leg_segment_clean(I,roi_img);
        
        bkg_sub = (ref_img - I);        
        I(:,:,2) = bkg_sub;
        
        I = I / 255;
        
        X{i_img,1} = I(:,:,1);
        X{i_img,2} = I(:,:,2);
        
        roi_img = imread([output_dir 'roi_' num2str(imgs_sec(i_img)) '.png']);
        roi_images{i_img} = padarray(roi_img,[20 20],'replicate'); 
    end
    
    score_images = batch_evaluate_boost_images(X,params,weak_learners,roi_images);
    
     for i_img = 1 : length(imgs_sec)
        I = imread([data_dir img_list(imgs_sec(i_img)).name]);        
        I = double(I)/255;
        I = padarray(I,[20 20],'replicate');       
        
        score_img = score_images{i_img};      
        est_leg = score_img > score_thres;
        %est_leg = bwdist(est_leg)<3 .*roi_images{i_img};
        %est_leg = est_leg & ~ neg_img{i_img} & roi_images{i_img};
        est_leg = filter_small_comp(est_leg,15);
    
        show_img_output = repmat(I,[1 1 3]);
        show_img_output(:,:,1) = show_img_output(:,:,1) + est_leg;
        imshow(show_img_output);
        
        imwrite(imcrop(show_img_output,[21 21 size(I,2)-41 size(I,1)-41]),[output_dir 'img_' num2str(imgs_sec(i_img)) '.png'],'png');
    end
end