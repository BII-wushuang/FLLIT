function [kernels,kernel_params,reg_tree] = remove_useless_filters(reg_tree,tmp_kernels,tmp_kernel_params)

% Remove the filters which have not been used in the weak-learner
%
%  authors: Carlos Becker, Roberto Rigamonti, CVLab EPFL
%  e-mail: name <dot> surname <at> epfl <dot> ch
%  web: http://cvlab.epfl.ch/
%  date: February 2014

var_idxs = [reg_tree.varIdx];
sub_wnd_list = unique(var_idxs(~[reg_tree.isLeaf]))+1;

kernels = cell(length(sub_wnd_list),1);
kernel_params = cell(length(sub_wnd_list),1);
for i_s = 1:length(sub_wnd_list)
    kernels{i_s} = tmp_kernels{sub_wnd_list(i_s)};
    kernel_params{i_s} = tmp_kernel_params{sub_wnd_list(i_s)};
    to_change = find([reg_tree.varIdx]==sub_wnd_list(i_s)-1);
    for i_c = 1:length(to_change)
        reg_tree(to_change(i_c)).varIdx = uint32(i_s-1);
    end
end

end
