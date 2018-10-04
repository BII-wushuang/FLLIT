function [raw_tips, raw_normtips] = findtips(imseg, imroi, locX, locY, centre_of_mass, theta)

imseg = imfill(imseg,'holes');
img = bwmorph(imseg,'thin','Inf');
imgnospurs = bwmorph(img,'spur',3);
% count the neighbours of the skeletons
neighbour_count=imfilter(uint8(imgnospurs),ones(3));
bw_junc = neighbour_count>3 & img;
bw_body = img - bw_junc;
bw_ends = bwmorph(bw_body,'endpoints');

connComp = bwconncomp(bw_body);

total_tips = 0;
for cc = 1:connComp.NumObjects
    clear endpoints;
    clear distmin;
    clear tipidx;
    
    legpixels = connComp.PixelIdxList{cc};

    imleg = zeros(size(imseg));
    imleg(legpixels) = 1;
    
    [row col] = find((imleg>0).*bw_ends);
    %[row col] = find((imleg>0));
    
    endpoints = [col row];
    
    npoints = length(endpoints(:,1));
    
    for j = 1:npoints
        [distmin(j),~] = min(pdist2([locX locY],endpoints(j,:),'euclidean'));
    end
    
    ntips = ceil(npoints/2);
    [~,tipidx] = sort(distmin);
    
    for i  =  1 : 1
        raw_tips(total_tips+i,:) = endpoints(tipidx(end+1-i),:);
    end
    
    total_tips = total_tips + ntips;
    
end

raw_normtips(:,1) = cosd(theta)*(raw_tips(:,1) - centre_of_mass(1)) + sind(theta)*(raw_tips(:,2) - centre_of_mass(2));
raw_normtips(:,2) = -sind(theta)*(raw_tips(:,1) - centre_of_mass(1)) + cosd(theta)*(raw_tips(:,2) - centre_of_mass(2));

end
