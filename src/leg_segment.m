function  [pos_img,neg_img_body,neg_img_bkg] = leg_segment(I,ref_img,sil_thres)
%pos_img: leg candidates on balck background
%neg_img_body: body of fruit fly withous legs
%neg_img_bkg: shade of fruit fly on white background

%ref_img is the background
if(nargin < 3)
    sil_thres = 0.1;
end

I = double(I);

fly_silhouette = max(ref_img - I,0) ./ ref_img > sil_thres;
%background segmentation, if the intensity ratio is greater than sil_thres,
%then keep it
fly_silhouette = bwareaopen(fly_silhouette, 200);
%return 1 if the pixel is above the threshold (500)
%fly_silhouette = imfill(fly_silhouette, 'holes');
neg_img_bkg = ones(size(I)) - fly_silhouette;
%return 1 on pixels with intensity less than threshold

cl_fly = bwmorph(fly_silhouette,'skel',Inf); %return skeleton structure of the image
edge_fly = bwmorph(fly_silhouette,'remove'); %return the outline of the image
%% consider this
%canny_edge_fly = edge(fly_silhouette, 'canny', 0.35)
%leg_candidate_1 = edge_fly .* (bwdist(cl_fly) < 2);
%leg_candidate_2 = canny_edge_fly .* (bwdist(cl_fly)<3);
%leg_candidate_1 = bwareaopen(leg_candidate_1,30);
%leg_candidate_2 = bwareaopen(leg_candidate_2,20);
%leg_candidate = leg_candidate_1 + leg_candidate_2;
%replace these line 29, 30 with the above codes
%%
leg_candidate = edge_fly .* (bwdist(cl_fly) < 3);
leg_candidate = bwareaopen(leg_candidate,25);

pos_img = (bwdist(leg_candidate) < 3) .* fly_silhouette;
pos_img = bwareaopen(pos_img,25);

neg_img_body = bwdist(gradientweight(double(fly_silhouette - pos_img)).*(fly_silhouette - pos_img))<6;
%gradientweight(double(fly_silhouette - pos_img)).*(fly_silhouette -
%pos_img) returns 1 near the body region (inc. body), thus this
%neg_img_body returns the area near the fly body

neg_img_body = bwareaopen(neg_img_body,200);

pos_img = (fly_silhouette-neg_img_body) > 0.5;
pos_img = bwareaopen(pos_img,10);
pos_img = (bwdist(pos_img) < 3) .* fly_silhouette;
pos_img = imfill(pos_img,'holes');

h = fspecial('average',[15 15]);
%create a filter with averaging effect of size 15*15 for smoothing
neg_img_body = imfilter(neg_img_body,h);
neg_img_body = neg_img_body > 0.7;
%neg_img_body returns the body of fruit fly and its surrounding area
