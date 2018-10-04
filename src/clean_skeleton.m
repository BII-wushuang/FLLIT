% A function to clean skeletons.
% 
function [bw_body,bw_junction,img]=clean_skeleton(seg)
%thin process
img=bwmorph(seg,'skel','Inf');
%img=bwmorph(img,'spur');
% count the neighbours of the skeletons
neighbour_count=imfilter(uint8(img),ones(3));
bw_body=neighbour_count<=3 & img;
bw_junction=neighbour_count>3 & img;
bw_ends=neighbour_count<=2 & img;
%find the terminal segments
bw_ends=imreconstruct(bw_ends,bw_body);

CC = bwconncomp(bw_ends);
numPixels = cellfun(@numel,CC.PixelIdxList);
[idx]=find(numPixels<3);

while ~isempty(idx)
    spur_length=2;
    % Remove the terminal segments if they are too short
    bw_body(bw_ends & ~bwareaopen(bw_ends, spur_length)) = false;
    img=(bw_body+bw_junction)>0;
    % Thin the binary image
    img = bwmorph(img, 'skel', Inf);
    %img=bwmorph(img,'spur');
    neighbour_count = imfilter(uint8(img), ones(3));
    bw_body = neighbour_count <=3 & img;
    bw_junction = neighbour_count >3 & img;
    bw_ends = neighbour_count <=2 & img;
    % Find the terminal segments - i.e. those containing end points
    bw_ends = imreconstruct(bw_ends, bw_body);
    
    CC = bwconncomp(bw_ends);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [idx]=find(numPixels<2);
end

img=(bw_body+bw_junction)>0;
img_no_soma=img;

%new_img_thin_wo_nu=bwmorph(new_img_thin_wo_nu,'spur');
neighbour_count = imfilter(uint8(img_no_soma), ones(3));
bw_body = neighbour_count <=3 & img_no_soma;
%figure,imshow(bw_body)
bw_junction = neighbour_count >3 & img_no_soma;
% This is for removing small segments in between two junctions.
bw_ends = neighbour_count <=2 & img_no_soma;
% Find the terminal segments - i.e. those containing end points
bw_ends = imreconstruct(bw_ends, bw_body);
bw_non_terminal=bw_body.*~(bw_ends); % non terminals

CC = bwconncomp(bw_non_terminal);
numPixels = cellfun(@numel,CC.PixelIdxList);
[idx]=find(numPixels<2);

for jj=1:length(idx)
    bw_junction(CC.PixelIdxList{idx(jj)})=true;
    bw_body(CC.PixelIdxList{idx(jj)})=false;
end

end