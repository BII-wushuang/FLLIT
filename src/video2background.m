function [thres, mean_I,std_I,ref_img,mask_imgs_fly] = video2background(data_dir, sub_dir) 

dir_result = './Results/Classifiers';
pos_bs = strfind(sub_dir, '/');

if (isempty(pos_bs))
    pos_bs = length(sub_dir);
end

if(~exist([dir_result sub_dir(1:pos_bs(end))]))
    mkdir([dir_result sub_dir(1:pos_bs(end))]);
end

m_img_fn = [dir_result  sub_dir '_Mean.mat'];
std_img_fn = [dir_result  sub_dir '_Std.mat'];

if (~isempty(dir([data_dir '*.tif'])))
    fn_list = dir([data_dir '*.tif']);
    %fly_sz = [100 100];
else
    fn_list = dir([data_dir '*.bmp']);
    %fly_sz = [250 250];
end

n_imgs = length(fn_list);

I = imread([data_dir fn_list(1).name]);

mean_I = zeros(size(I)+40);

if(exist(m_img_fn))
    load(m_img_fn);
    load(std_img_fn);
else
    
    for i_img = 1 : length(fn_list)
        I = imread([data_dir fn_list(i_img).name]);
        I = double(I);
        I = padarray(I,[20 20],'replicate');
        mean_I = mean_I + I;
    end
    
    mean_I = mean_I / n_imgs;
    std_I = zeros(size(mean_I));
    
    for i_img = 1 : length(fn_list)
        I = imread([data_dir fn_list(i_img).name]);
        I = double(I);
        I = padarray(I,[20 20],'replicate');
        std_I = abs(I - mean_I) + std_I;
    end
    
    std_I = std_I / n_imgs;
    save(m_img_fn,'mean_I');
    save(std_img_fn,'std_I');
end

% obtain the background image
mask_imgs_fly = zeros(size(mean_I));
bkg_imgs = zeros(size(mean_I));

thres = 0.1;

sample_ratio = 25;

for i_img = 1 : floor(length(fn_list) / sample_ratio)
    I = imread([data_dir fn_list(i_img * sample_ratio).name]);
    I = double(I);
    I = padarray(I,[20 20],'replicate');
    tmp_mask =  (max(mean_I - I,0) ./ mean_I) < thres;
    
    mask_imgs_fly = mask_imgs_fly + tmp_mask;
    bkg_imgs = bkg_imgs + I .* tmp_mask;
end

ref_img = bkg_imgs ./ mask_imgs_fly;