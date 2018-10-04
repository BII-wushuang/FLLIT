%Classifier for leg segmentation
function Segmentation_Base (data_dir)

if(nargin<1)
    data_dir = uigetdir('./Data');
end

pos_bs = strfind(data_dir,'Data');
sub_dir = data_dir(pos_bs(end)+length('Data'):length(data_dir));
data_dir = ['./Data' sub_dir '/'];
output_dir = ['./Results/SegmentedImages' sub_dir '/'];
if(~exist(output_dir))
    mkdir(output_dir);
end

[thres, ~,~,ref_img,~] = video2background(data_dir, sub_dir);

if(~isempty(dir([data_dir '*.tif'])))
    fn_list = dir([data_dir '*.tif']);
else
    fn_list = dir([data_dir '*.bmp']);
end

for i_img = 1 : length(fn_list)
    I = imread([data_dir fn_list(i_img).name]);
    I = double(I);
    I = padarray(I,[20 20],'replicate');

    roi_img = (max(ref_img - I(:,:,1),0) ./ ref_img) > thres;

    [pos_img,~,~] = leg_segment(I,ref_img,thres);
    %pos_img = (filter_small_comp(pos_img,15) - imroi)>0.1;

    show_img_output = repmat(double(I)/255,[1 1 3]);
    show_img_output(:,:,1) = show_img_output(:,:,1) + pos_img;

    imwrite(imcrop(show_img_output,[21 21 size(I,2)-41 size(I,1)-41]),[output_dir 'img_' num2str(i_img) '.png'],'png');
    imwrite(imcrop(roi_img,[21 21 size(I,2)-41 size(I,1)-41]),[output_dir 'roi_' num2str(i_img) '.png'],'png');
end

end