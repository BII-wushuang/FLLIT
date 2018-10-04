function [good_sampling_idx,tot_idx] = get_sampling_data(params,data_train)

% Find good positions for sampling
%
%  authors: Carlos Becker, Roberto Rigamonti, CVLab EPFL
%  e-mail: name <dot> surname <at> epfl <dot> ch
%  web: http://cvlab.epfl.ch/
%  date: February 2014

imgs_no = size(data_train.imgs.X,1);


good_sampling_idx = cell(imgs_no,1);
selected_regions = cell(imgs_no,1);
for i_gs = 1:imgs_no
    selected_regions{i_gs} = data_train.masks.X{i_gs};
end

% Get rid of image borders
for i_gs = 1:imgs_no
    border_img = zeros(size(selected_regions{i_gs}));
    border_img(1:params.border_skip_size,:) = 1;
    border_img(:,1:params.border_skip_size) = 1;
    border_img((end-params.border_skip_size):end,:) = 1;
    border_img(:,(end-params.border_skip_size):end) = 1;
    selected_regions{i_gs} = selected_regions{i_gs} & ~border_img;
end

tot_idx.pos = 0;
tot_idx.neg = 0;
% Take as positive indices the indices in the GT, the rest as negatives
for i_gs = 1:imgs_no
    good_sampling_idx{i_gs}.pos = find(data_train.gts.X{i_gs} & selected_regions{i_gs});
    good_sampling_idx{i_gs}.neg = find(~data_train.gts.X{i_gs} & selected_regions{i_gs});
    tot_idx.pos = tot_idx.pos+length(good_sampling_idx{i_gs}.pos);
    tot_idx.neg = tot_idx.neg+length(good_sampling_idx{i_gs}.neg);
end

end
