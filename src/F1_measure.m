function [F,PPV,TPR,Pre,Rec,Spec,MCC,thres,Accu,Jaccard_Index,AP] = F1_measure(E,gt_I,mask_I, weight_image,thres)

% measure the segmenation performance using F measure

% input: E         Estimated segmentation
%        gt_I      grountruth segmentation
%        mask_I    only pixels in the mask_I are accounted
%        thresh    assigned threshold, if the threshold is not given,
%        output the result with highest F score

% output: F       Fmeasure
%         PPV     Precision used for PR-curve
%         TPR     Recall used for PR-curve
%         Pre     Precision
%         Rec     Recall
%         Spec    Specificity
%         MCC     Mattews correlation coefficient
%         Accu    Accuracy
%         Jaccard_Index  Jaccard Index also kown as VOC score
%         AP      Average precision
%         thres   the selected threhold


if(nargin < 5)

    thres = 0.01: (1 / 256) : 1;
    
end

if(nargin < 4)

    weight_image = ones(size(E));
    
end


if(nargin < 3)

    mask_I = ones(size(E));
    
end


E = E .* mask_I;

gt_I = gt_I .* mask_I;


for i_t = 1 : length(thres)
    
    E1 = E > thres(i_t);  
    
    E1 = E1 .* mask_I;
        
    PPV(i_t) = sum(sum(E1 .* gt_I .* weight_image)) / max(1,sum(E1(:) .* weight_image(:)));
    
    TPR(i_t) =  sum(sum(E1 .* gt_I .* weight_image)) / max(1,sum(gt_I(:) .* weight_image(:)));
    
end

%[F,thre_idx] = max(2* PPV .* TPR ./ (PPV + TPR));
[F,thre_idx] = max(1.25* PPV .* TPR ./ (0.25*PPV + TPR));

% thres_given = thres(thre_idx);

E_I = E > thres(thre_idx);

E_I = E_I .* mask_I;

TP = sum(gt_I(:) .* (E_I(:) .* weight_image(:)));
TN = sum(~gt_I(:) .* (~E_I(:) .* weight_image(:)));
FP = sum(~gt_I(:) .* (E_I(:) .* weight_image(:)));
FN = sum(gt_I(:) .* (~E_I(:) .* weight_image(:)));
N = sum(mask_I(:));

Pre = PPV(thre_idx);
Rec = TPR(thre_idx);
Spec = TN  / (TN + FP);
MCC = (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN));

% S = (TP+FN) / N;
% P = (TP+FP)/N;
% MCC1 = (TP / N - S *P) / sqrt(P*S*(1-S)*(1-P));

Accu = sum(gt_I(:) .* E_I(:) .* weight_image(:) + (~gt_I(:) .* (~E_I(:) .* weight_image(:)))) / (size(gt_I,1) * size(gt_I,2)) / sum(weight_image(:));

Jaccard_Index = sum(gt_I(:) .* E_I(:) .* weight_image(:)) / (sum(gt_I(:) .* weight_image(:)) + sum(E_I(:) .* weight_image(:)) - sum(gt_I(:) .* E_I(:) .* weight_image(:)));

thres = thres(thre_idx);

AP = mean(PPV);
