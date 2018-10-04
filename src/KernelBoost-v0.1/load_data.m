function [params,data] = load_data(params)

% Load program's data
%
%  authors: Carlos Becker, Roberto Rigamonti, CVLab EPFL
%  e-mail: name <dot> surname <at> epfl <dot> ch
%  web: http://cvlab.epfl.ch/
%  date: February 2014

for i_op = 1:params.ops_no
    op = params.ops{i_op};
    
    %% Load GTs and masks
    data.(op).gts.X = load_binary_channel(params,op,'gts');
    data.(op).masks.X = load_binary_channel(params,op,'masks');
    
    %% Load the remaining channels
    for i_ch = 1:params.ch_no
        ch = params.ch_list{i_ch};
        
        switch ch
            case 'imgs'
                data.(op).(ch).X = load_plain_data(params,op,ch);
                data.(op).(ch).idxs = 1;
                data.(op).(ch).sub_ch_no = 1;
            otherwise
                error('Loading function for channel %s not implemented yet!',ch);
        end
    end
end

end
