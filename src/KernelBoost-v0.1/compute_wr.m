function [wrIdxs,wrResponses,wrWeights] = compute_wr(params,T1_idx,W,R,computeIndivLoss,compute2ndDeriv,labels,current_response)

% Compute weights and responses
%
%  authors: Carlos Becker, Roberto Rigamonti, CVLab EPFL
%  e-mail: name <dot> surname <at> epfl <dot> ch
%  web: http://cvlab.epfl.ch/
%  date: February 2014

if (params.use_uniform_random_sampling&&params.use_2nd_deriv)
    error('Conflicting flags: uniform and 2nd derivative flags activated!');
end

if (params.use_qws)
    samplingWeights = computeIndivLoss(current_response);
elseif (params.use_2nd_deriv)
    samplingWeights = compute2ndDeriv(current_response);
else
    samplingWeights = W;
end

if (params.use_uniform_random_sampling)
    samplingWeights(:) = 1/length(samplingWeights);
end

if (params.use_qws)
    newAll = params.pos_to_sample_no+params.neg_to_sample_no;
    wrIdxs = mexQws(samplingWeights(T1_idx),uint32(newAll));
    wrIdxs = T1_idx(unique(wrIdxs));
    wrResponses = R(wrIdxs);
    wrWeights = W(wrIdxs);
else
    posIdxs = T1_idx(labels(T1_idx)>0)';
    negIdxs = T1_idx(labels(T1_idx)<0)';

    % sample half at random
    newNNeg = round(length(negIdxs)/2);
    newNPos = round(length(posIdxs)/2);
    
    newPos = randsample(posIdxs,newNPos);
    newNeg = randsample(negIdxs,newNNeg);

    wrIdxs = [newPos;newNeg];
    wrLabels = labels(wrIdxs);
    wrResponses = R(wrIdxs);

    % normalize weights btw pos/neg
    negWeightsMul = sum(W(negIdxs))/sum(W(newNeg));
    posWeightsMul = sum(W(posIdxs))/sum(W(newPos));
    wrWeights = W(wrIdxs);
    wrWeights(wrLabels>0) = wrWeights(wrLabels>0)*posWeightsMul;
    wrWeights(wrLabels<0) = wrWeights(wrLabels<0)*negWeightsMul;
end
% Normalize weights
wrWeights = wrWeights/sum(wrWeights);

end
