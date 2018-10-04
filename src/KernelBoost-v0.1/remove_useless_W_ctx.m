function [kernels,kernel_params,reg_tree,sub_wnd_kern_list,ctx_list] = remove_useless_W_ctx(reg_tree,tmp_W,tmp_W_info)

% Remove the W which have not been used in the weak-learner
%
%  authors: Carlos Becker, Roberto Rigamonti, CVLab EPFL
%  e-mail: name <dot> surname <at> epfl <dot> ch
%  web: http://cvlab.epfl.ch/
%  date: February 2014

var_idxs = [reg_tree.varIdx];
sub_wnd_list = unique(var_idxs(~[reg_tree.isLeaf]))+1;

% kernels = cell(length(sub_wnd_list),1);
% kernel_params = cell(length(sub_wnd_list),1);


k_idx = 1;

for ik = 1 : length(tmp_W)
   
    for i_w = 1 : size(tmp_W{ik},2)
    
        tmp_kernels{k_idx} = tmp_W{ik}(:,i_w);
        
        tmp_kernel_params{k_idx} = tmp_W_info(ik,:);
        
        k_idx = k_idx + 1;
    
    end
    
end


n_kernels = length(tmp_kernels);

sub_wnd_kern_list = sub_wnd_list;

sub_wnd_kern_list(sub_wnd_list > n_kernels) = [];

n_var_nk = sum(sub_wnd_list > n_kernels);

n_var_k = length(sub_wnd_list) - n_var_nk;

kernels = cell(n_var_k,1);

kernel_params = cell(n_var_k,1);

if(n_var_nk > 0)
    
    ctx_list = sub_wnd_list(n_var_k + 1 : end) - n_kernels;
    
else
    
    ctx_list = [];
    
end


for i_s = 1:length(sub_wnd_list)
    
    if(sub_wnd_list(i_s) < (n_kernels + 1))
    
        kernels{i_s} = tmp_kernels{sub_wnd_list(i_s)};
    
        kernel_params{i_s} = tmp_kernel_params{sub_wnd_list(i_s)};
        
    end
    
    to_change = find([reg_tree.varIdx]==sub_wnd_list(i_s)-1);
    
    for i_c = 1:length(to_change)
        
        reg_tree(to_change(i_c)).varIdx = uint32(i_s-1);
        
    end
end

end
