function [reg_tree,sub_wnd_list] = remove_ftrs(reg_tree)

% Remove the filters which have not been used in the weak-learner
%
%  authors: Carlos Becker, Roberto Rigamonti, CVLab EPFL
%  e-mail: name <dot> surname <at> epfl <dot> ch
%  web: http://cvlab.epfl.ch/
%  date: February 2014

var_idxs = [reg_tree.varIdx];
sub_wnd_list = unique(var_idxs(~[reg_tree.isLeaf]))+1;

for i_s = 1:length(sub_wnd_list)
    
    to_change = find([reg_tree.varIdx]==sub_wnd_list(i_s)-1);
    
    for i_c = 1:length(to_change)
    
        reg_tree(to_change(i_c)).varIdx = uint32(i_s-1);
    
    end
    
end

end
