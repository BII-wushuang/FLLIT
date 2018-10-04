function [kernels,kernel_params,features] = merge_features_kernels(tmp_kernels,tmp_kernel_params,tmp_features)

% Merge kernels, params and features from the different
% channels/subchannels
%
%  authors: Carlos Becker, Roberto Rigamonti, CVLab EPFL
%  e-mail: name <dot> surname <at> epfl <dot> ch
%  web: http://cvlab.epfl.ch/
%  date: February 2014

k_count = 0;
for i_ch = 1:length(tmp_kernel_params)
    if (length(tmp_kernel_params{i_ch})>0)
        k_count = k_count+length(tmp_kernel_params{i_ch})*length(tmp_kernel_params{i_ch}{1});
    end
end

kernels = cell(k_count,1);
kernel_params = cell(k_count,1);
features = zeros(size(tmp_features{1}{1},1),k_count);

i_count = 1;
for i_ch = 1:length(tmp_kernel_params)
    for i_s = 1:length(tmp_kernel_params{i_ch})
        kernels(i_count:i_count+length(tmp_kernel_params{i_ch}{i_s})-1) = tmp_kernels{i_ch}{i_s}(:);
        kernel_params(i_count:i_count+length(tmp_kernel_params{i_ch}{i_s})-1) = tmp_kernel_params{i_ch}{i_s}(:);
        features(:,i_count:i_count+length(tmp_kernel_params{i_ch}{i_s})-1) = tmp_features{i_ch}{i_s};
        i_count = i_count+length(tmp_kernel_params{i_ch}{i_s});
    end
end

end
