function [kernels,kernel_params,reg_tree,sub_wnd_list] = remove_useless_W(reg_tree,tmp_W,tmp_W_info)

% Remove the W which have not been used in the weak-learner
%
%  authors: Carlos Becker, Roberto Rigamonti, CVLab EPFL
%  e-mail: name <dot> surname <at> epfl <dot> ch
%  web: http://cvlab.epfl.ch/
%  date: February 2014

var_idxs = [reg_tree.varIdx];
sub_wnd_list = unique(var_idxs(~[reg_tree.isLeaf]))+1;

kernels = cell(length(sub_wnd_list),1);
kernel_params = cell(length(sub_wnd_list),1);


k_idx = 1;

for ik = 1 : length(tmp_W)
   
    for i_w = 1 : size(tmp_W{ik},2)
    
        tmp_kernels{k_idx} = tmp_W{ik}(:,i_w);
        
        tmp_kernel_params{k_idx} = tmp_W_info(ik,:);
        
        k_idx = k_idx + 1;
    
    end
    
end





for i_s = 1:length(sub_wnd_list)
    kernels{i_s} = tmp_kernels{sub_wnd_list(i_s)};
    kernel_params{i_s} = tmp_kernel_params{sub_wnd_list(i_s)};
    to_change = find([reg_tree.varIdx]==sub_wnd_list(i_s)-1);
    for i_c = 1:length(to_change)
        reg_tree(to_change(i_c)).varIdx = uint32(i_s-1);
    end
end

end
