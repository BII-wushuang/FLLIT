function [normalized] = normalize_min_max_mask(in_matrix,mask)

scan_area = in_matrix(mask>0);
min_val = min(scan_area(:));
max_val = max(scan_area(:));
if (abs(max_val-min_val)>1e-5)
    normalized = (in_matrix-min_val)/(max_val-min_val);
else
    normalized = zeros(size(in_matrix));
end

normalized(mask==0) = 0;

end
