function Image = connectify(img, imroi, imseg)
%CONNECTIFY Find the disconnected components and group them leg-wise. 
% Convert the pixel level segmentation result to an object level.

se1 = strel('disk',5);
se2 = strel('disk',0);

% Morphological operations to estimate the body of the fly
imbody = imerode(imroi,se1);
imbody = imdilate(imbody,se1);
% %figure, imshow(imfuse(imroi,imbody));

if (0)
    imsk_old = bwmorph(img,'skel',Inf)*255;
    imsk_new = bwmorph(skeleton(img) > 10, 'skel', Inf) * 255;
    montage_image(:,:,1,1) = imsk_old;
    montage_image(:,:,1,2) = imsk_new;
    figure, montage(montage_image);
end

connComp = bwconncomp(img);
nRegions = connComp.NumObjects;
% Get information on the segmentation's centroid, medial axis with its 
% end points and junctions. This will then be used to
% identify whether regions must be connected or not!
for cc = 1:nRegions
    legpixels = connComp.PixelIdxList{cc};
    [I,J] = ind2sub(size(img),legpixels);
    connComp.centroid{cc} = [mean(I) mean(J)];

    imleg = zeros(size(img));
    imleg(legpixels) = 1;
    [skr,rad] = skeleton(imleg);
    imsk = bwmorph(skr > 10,'skel',Inf);
    % Skeleton analysis: Get the end points.
    [dmap,exy,jxy] = anaskel(imsk);

    connComp.exy{cc} = exy;
    connComp.jxy{cc} = jxy;
    connComp.medialAxis{cc} = find(imsk > 0);
end

distThresh = 10; %17.5

PT1 = [];
PT2 = [];
for i = 1:nRegions
    for j = i+1:nRegions
        exyI = connComp.exy{i};
        exyJ = connComp.exy{j};
        [p1_,p2_] = checkToConnect(exyI,exyJ,distThresh);
        PT1 = [PT1 p1_];
        PT2 = [PT2 p2_];
    end
end

Image = img;

npoints = size(PT1,2);
for i = 1:npoints    
    x = [PT1(1,i) PT2(1,i)];
    y = [PT1(2,i) PT2(2,i)];
    
    Image = func_DrawLine(Image,y(1),x(1),y(2),x(2),1);
    %X = PT1(1,i):PT2(1,i);
    %Y = round(interp1(x,y,X));    
    %ind = sub2ind(size(Image),Y,X);
    %Image(ind) = 255;
end

end

function [pt1,pt2] = checkToConnect(exyI,exyJ,distThresh)

pt1 = [];
pt2 = [];

nendpoints_i = length(exyI(1,:));
for k = 1:nendpoints_i
    A = bsxfun(@minus,exyJ,exyI(:,k));
    distances = sqrt(sum(A.^2));
    
    idx = find(distances < distThresh);
    for i = 1:length(idx)
        % exyJ(:,idx(i)) and exyI(:,k) are close endpoints of some two medial
        % axis. Now must enforce directionality constraint
        probableLeg1 = exyI(:,k) - exyI(:,~(k-1)+1);
        probableLeg1 = probableLeg1./norm(probableLeg1);
        
        %probableLeg2 = exyJ(:,~(idx(i)-1)+1) - exyJ(:,idx(i));
        probableJoin = exyJ(:,idx(i)) - exyI(:,k);
        probableJoin = - (probableJoin./norm(probableJoin));
        
        theta = acosd(probableLeg1' * probableJoin);
        
        if (theta > 120) %|| (theta < 30)
            % If this also holds, then mark this point
            pt1 = [pt1 exyI(:,k)];
            pt2 = [pt2 exyJ(:,idx(i))];
        end
    end
end

end