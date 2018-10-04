function [F,Pre,Rec,Accu,Jaccard_Index,AP,thres] = seg_measure_morph(E,gt_I,mask_I,thres)

% measure the segmenation performance, here a morphology is used to improve
% the performance

% input: E         Estimated segmentation
%        gt_I      grountruth segmentation
%        mask_I    only pixels in the mask_I are accounted
%        thresh    assigned threshold, if the threshold is not given,
%        output the result with highest F score

% output: F       Fmeasure
%         Pre     Precision
%         Rec     Recall
%         Accu    Accuracy
%         Jaccard_Index  Jaccard Index also kown as VOC score
%         AP      Average precision
%         thres   the selected threhold



if(nargin < 4)

    thres = 0.01: 0.01 : 1;
    
end
    

    gt_I = gt_I .* mask_I;
    
for i_t = 1 : length(thres)
    
    E1 = E > thres(i_t);
    
    CC = bwconncomp(E1);
    
    numPixels = cellfun(@numel,CC.PixelIdxList);
    
    idx_n = find(numPixels < 200);
    
    for i_idx = 1 : length(idx_n)
        
        E1(CC.PixelIdxList{idx_n(i_idx)}) = 0;
        
    end
    
    
    E1 = E1 .* mask_I;
    

    
    PPV(i_t) = sum(sum(E1 .* gt_I)) / max(1,sum(E1(:)));
    
    TPR(i_t) =  sum(sum(E1 .* gt_I)) / max(1,sum(gt_I(:)));
    
    
end

[F,thre_idx] = max(PPV .* TPR ./ (PPV / 2 + TPR / 2));

% thres_given = thres(thre_idx);

Pre = PPV(thre_idx);

Rec = TPR(thre_idx);

E_I = E > thres(thre_idx);

CC = bwconncomp(E_I);

numPixels = cellfun(@numel,CC.PixelIdxList);

idx_n = find(numPixels < 200);

for i_idx = 1 : length(idx_n)
    
    E_I(CC.PixelIdxList{idx_n(i_idx)}) = 0;
    
end


E_I = E_I .* mask_I;



Accu = sum(gt_I(:) .* E_I(:) + (~gt_I(:) .* (~E_I(:)))) / (size(gt_I,1) * size(gt_I,2));

Jaccard_Index = sum(gt_I(:) .* E_I(:)) / (sum(gt_I(:)) + sum(E_I(:)) - sum(gt_I(:) .* E_I(:)));

thres = thres(thre_idx);

AP = mean(PPV);

    
    
% else
%     
%     
%     
% end