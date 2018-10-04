function [compute_wi,compute_ri,compute_loss,compute_indiv_loss,compute_2nd_deriv,mex_loss_type] = select_fncts(params,labels)

% Select the proper functions according to the chosen loss
%
%  authors: Carlos Becker, Roberto Rigamonti, CVLab EPFL
%  e-mail: name <dot> surname <at> epfl <dot> ch
%  web: http://cvlab.epfl.ch/
%  date: February 2014

switch params.loss_type
    case 'exp'
        compute_wi = @(resp) exp(-labels.*resp);
        compute_ri = @(resp) labels;
        compute_loss = @(resp) sum(exp(-labels.*resp))/length(labels);
        compute_indiv_loss = @(resp) exp(-labels.*resp)/length(labels);
        compute_2nd_deriv = @(resp) compute_wi(resp);
        mex_loss_type = 'exploss';
    case 'expWithMAdaBoost'
        compute_wi = @(resp) max(exp(-labels.*resp),1.0);
        compute_ri = @(resp) labels;
        compute_loss = @(resp) sum(max(exp(-labels.*resp),1.0))/length(labels);
        compute_indiv_loss = @(resp) compute_wi(resp);
        compute_2nd_deriv = @(resp) compute_wi(resp);
        mex_loss_type = 'exploss';
    case 'log'
        % this is not SQB, but gradient boosting, for numerical stability
        compute_wi = @(resp) ones(size(resp,1),1);
        compute_ri = @(resp) (2*exp(-2*labels.*resp).*labels)./(1+exp(-2*labels.*resp));
        compute_loss = @(resp) sum(log(eps+1+exp(-2*labels.*resp)))/length(labels);
        compute_indiv_loss = @(resp) log(eps+1+exp(-2*labels.*resp))/length(labels);
        compute_2nd_deriv = @(resp) 4*exp(-2*labels.*resp)./(1+exp(-2*labels.*resp)).^2;
        mex_loss_type = 'logloss';
    case 'squared'
        compute_wi = @(resp) 2*ones(size(resp,1),1);
        compute_ri = @(resp) labels-resp;
        compute_loss = @(resp) sum((resp-labels).^2)/length(labels);
        compute_indiv_loss = @(resp) (resp-labels).^2/length(labels);
        compute_2nd_deriv = [];
        mex_loss_type = 'squaredloss';
    case 'RandomForest'
        % this is random forest
        compute_wi = @(resp) ones(size(resp,1),1);
        compute_ri = @(resp) labels;
        compute_loss = @(resp) sum(exp(-labels.*resp))/length(labels);
        compute_indiv_loss = @(resp) exp(-labels.*resp)/length(labels);
        compute_2nd_deriv = @(resp) compute_wi(resp);
        mex_loss_type = 'exploss';
    otherwise
        error('Incorrect loss type specified!');
end

end
