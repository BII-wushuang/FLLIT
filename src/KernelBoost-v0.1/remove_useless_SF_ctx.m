function [reg_tree,sub_wnd_kern_list,ctx_list] = remove_useless_SF_ctx(reg_tree,n_kernels)

% Remove the features which have not been used in the weak-learner
%
%  authors: Carlos Becker, Roberto Rigamonti, CVLab EPFL
%  e-mail: name <dot> surname <at> epfl <dot> ch
%  web: http://cvlab.epfl.ch/
%  date: February 2014

var_idxs = [reg_tree.varIdx];

sub_wnd_list = unique(var_idxs(~[reg_tree.isLeaf]))+1;


sub_wnd_kern_list = sub_wnd_list;

sub_wnd_kern_list(sub_wnd_list > n_kernels) = [];

n_var_nk = sum(sub_wnd_list > n_kernels);

n_var_k = length(sub_wnd_list) - n_var_nk;

if(n_var_nk > 0)
    
    ctx_list = sub_wnd_list(n_var_k + 1 : end) - n_kernels;
    
else
    
    ctx_list = [];
    
end


for i_s = 1:length(sub_wnd_list)
    
    to_change = find([reg_tree.varIdx]==sub_wnd_list(i_s)-1);
    
    for i_c = 1:length(to_change)
        
        reg_tree(to_change(i_c)).varIdx = uint32(i_s-1);
        
    end
end

end
