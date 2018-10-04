function [params] = setup_config_L_v2(dataset_name)

% Choose the configuration file according to the selected dataset
%
%  authors: Carlos Becker, Roberto Rigamonti, CVLab EPFL
%  e-mail: name <dot> surname <at> epfl <dot> ch
%  web: http://cvlab.epfl.ch/
%  date: February 2014

switch dataset_name
    case {'DRIVE'}
        params = setup_config_DRIVE_L(dataset_name);
    case {'STARE'}
        params = setup_config_STARE(dataset_name);
            case {'HRF'}
        params = setup_config_HRF(dataset_name);
        
    case {'Neuron'}
        
         params = setup_config_DRIVE_L(dataset_name);
        
        
    otherwise
%         error('Dataset not recognized (%s)',dataset_name);

params = setup_config_DRIVE_L(dataset_name);
end

end
