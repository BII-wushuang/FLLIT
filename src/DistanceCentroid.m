data_dirs = dir('./Data/Centroid/');
fileID = fopen(['DistanceCentroid.csv'],'w');
fprintf(fileID,'%s \t', 'Image', 'Uncorrected Length', 'Distance to Centroid', 'Ratio');
fprintf(fileID,'\n');

for i = 3:size(data_dirs,1)
    data_dir = ['./Data/Centroid/' data_dirs(i).name];
    pos_bs = strfind(data_dir,'Data');
    sub_dir = data_dir(pos_bs(end)+length('Data'):length(data_dir));
    data_dir = ['./Data' sub_dir '/'];
    
    output_dir = ['./Results/SegmentedImages' sub_dir '/'];
    if(~exist(output_dir))
        mkdir(output_dir);
    end
    
    fn_list = dir([data_dir '*.tif']);
%     I = imread([data_dir fn_list(1).name]);
%     I = rgb2gray(I(:,:,1:3));
%     [thres,~,~,ref_img,~] = video2background(data_dir, sub_dir);
%     imshow(imcrop(ref_img,[21 21 size(I,2)-1 size(I,1)-1]),[]);
%     
%     for i_img = 1 : length(fn_list)
%         I = imread([data_dir fn_list(i_img).name]);
%         I = rgb2gray(I(:,:,1:3));
%         I = double(I);
%         I = padarray(I,[20 20],'replicate');
% 
%         roi_img = (max(ref_img - I(:,:,1),0) ./ ref_img) > thres;
%         roi_img = bwareaopen(roi_img, 1000);
%         imwrite(imcrop(roi_img,[21 21 size(I,2)-41 size(I,1)-41]),[output_dir 'roi_' num2str(i_img) '.png'],'png');
%     end
    
    se = strel('disk',4);
    for i_img = 1 : length(fn_list)
        I = imread([data_dir fn_list(i_img).name]);
        I = I(:,:,1:3);

        roi_img = imread([output_dir 'roi_' num2str(i_img) '.png']);
        Imroi = imerode(roi_img,se);
        [locY,locX] = find(Imroi > 0);
        [CoM,theta] = findCoM(locX,locY);
        
        locX_norm = locX - CoM(1);
        locY_norm = locY - CoM(2);
        points = [cosd(theta) sind(theta); -sind(theta) cosd(theta)]* [locX_norm'; locY_norm'];

        %Attempt to ensure that the head is in the positive
        %y-direction by comparing which side has more number of
        %pixel points
        if (length(find(points(2,:) > 0)) < length(find(points(2,:) < 0)))
            theta = mod(theta + 180,360);
        end
        
        img_norm = imtranslate(roi_img, [255 - CoM(1), 255 - CoM(2)]);
        img_norm = imrotate(img_norm, theta);
        img_norm = imcrop(img_norm, [size(img_norm,1)/2-150 size(img_norm,2)/2-150 300 300]);
        [Y,X] = find(img_norm);
        len = max(Y(find(X==150)))-min(Y(find(X==150)));
        
        [rows, columns] = find(I(:, :, 1) == 237);
        new_I = zeros(512,512);
        new_I(rows, columns) = 1;
        
        I_norm = imtranslate(new_I, [255 - CoM(1), 255 - CoM(2)]);
        I_norm = imrotate(I_norm, theta);
        I_norm = imcrop(I_norm, [size(I_norm,1)/2-150 size(I_norm,2)/2-150 300 300]);
        [rows, columns] = find(I_norm);
        dist_CoM = abs(150 - rows(1));
        
        fprintf(fileID, '%s \t', fn_list(i_img).name(1:end-4));
        fprintf(fileID, '%d \t', len, dist_CoM, dist_CoM/len);
        fprintf(fileID, '\n');
    end
end

fclose(fileID);
