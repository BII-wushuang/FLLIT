function BW = filter_small_comp(BW,t)

CC = bwconncomp(BW);

for ic = 1 : CC.NumObjects
    
    PIL = CC.PixelIdxList{ic};
    
    if(length(PIL) < t)
        
        BW(PIL) = 0;
        
    end
    
end
