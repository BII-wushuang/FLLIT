clear all; close all; clc;

data_dir = uigetdir('./Results/Tracking/');
load([data_dir '\CoM.mat']);
pos = strfind(data_dir,'Tracking');
sub_dir = data_dir(pos+length('Tracking'):length(data_dir));
pathRoiResults = ['./Results/SegmentedImages' sub_dir];

start_frame = find(CoM,1);
end_frame = length(dir([pathRoiResults '\*.png']))/2;

for i = start_frame: end_frame
    img = imread([pathRoiResults '\roi_' num2str(i) '.png']);
    img_norm = imtranslate(img, [255 - CoM(i,1), 255 - CoM(i,2)]);
    img_norm = imrotate(img_norm, CoM(i,3));
    img_norm = imcrop(img_norm, [size(img_norm,1)/2-150 size(img_norm,2)/2-150 300 300]);
    [Y,X] = find(img_norm);
    length(i) = max(Y(find(X==150)))-min(Y(find(X==150)));
end

mean(length(find(length>0)))
std(length(find(length>0)))