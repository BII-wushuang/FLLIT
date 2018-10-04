%%Centroid Measurements
folders = dir(['./Centroid']);

data_dirs = cell([length(folders)-2,1]);
for i = 1 : length(data_dirs)
    data_dirs{i} = ['./Centroid/' folders(i+2).name];
    if(~exist([data_dirs{i} '/roi']))
        mkdir([data_dirs{i} '/roi']);
        fn_list = dir([data_dirs{i} '/*.tif']);
        n_imgs = length(fn_list);
        I = imread([data_dirs{i} '/' fn_list(1).name]);
        mean_I = zeros(size(I));

        % background processing
        for i_img = 1 : length(fn_list)
            I = imread([data_dirs{i} '/' fn_list(i_img).name]);
            I = double(I);
            mean_I = mean_I + I(:,:,1);
        end

        mean_I = mean_I / n_imgs;
        std_I = zeros(size(mean_I));

        for i_img = 1 : length(fn_list)
            I = imread([data_dirs{i} '/'  fn_list(i_img).name]);
            I = double(I);
            std_I = abs(I(:,:,1) - mean_I) + std_I;
        end

        std_I = std_I / n_imgs;

        % obtain the background image
        mask_imgs_fly = zeros(size(mean_I));
        bkg_imgs = zeros(size(mean_I));

        thres = 0.1;

        sample_ratio = 25;
        
        for i_img = 1 : floor(length(fn_list) / sample_ratio)
            I = imread([data_dirs{i} '/'  fn_list(i_img * sample_ratio).name]);
            I = double(I);
            tmp_mask =  (max(mean_I - I(:,:,1),0) ./ mean_I) < thres;

            mask_imgs_fly = mask_imgs_fly + tmp_mask;
            bkg_imgs = bkg_imgs + I(:,:,1) .* tmp_mask;
        end

        ref_img = bkg_imgs ./ mask_imgs_fly;
        
        % background subtraction
        for i_img = 1 : length(fn_list)
            I = imread([data_dirs{i} '/'  fn_list(i_img).name]);
            I = double(I);
            fly_silhouette = max(ref_img - I(:,:,1),0) ./ ref_img > 0.2;
            fly_silhouette = bwareaopen(fly_silhouette,500);
            imwrite(fly_silhouette,[data_dirs{i} '/roi/roi_' num2str(i_img) '.png']);
        end
    end
    
    se = strel('disk',4);
    if(~exist([data_dirs{i} '/CoM.csv']))
         fn_list = dir([data_dirs{i} '/*.tif']);
         CoM = zeros([length(fn_list),3]);
         for i_img = 1 : length(fn_list)
            Imroi = imread([data_dirs{i} '/roi/roi_'  num2str(i_img) '.png']);
            Imroi = imerode(Imroi,se);
            [locY,locX] = find(Imroi > 0);
            [CoM(i_img,1:2) CoM(i_img,3)] = findCoM(locX,locY);
            if(i_img > 1)
                if (abs(CoM(i_img,3) - CoM(i_img-1,3)) > 90 && abs(CoM(i_img,3) - CoM(i_img-1,3)) < 270)
                    CoM(i_img,3) = mod(CoM(i_img,3) + 180,360);
                end
            end
         end
         csvwrite([data_dirs{i} '/CoM.csv'],CoM);
    end
end

for i = 1 : length(data_dirs)
    imgs = dir([data_dirs{i} '/New folder/*.tif']);
    CoM = csvread([data_dirs{i} '/CoM.csv']);
    for j = 1:size(imgs,1)
        imori = imread([data_dirs{i} '/New folder/' imgs(j).name]);
        imori = imori(:,:,1:3);
        [a,b] = find(imori(:,:,1)==237);
        i_img = str2num(imgs(j).name(end-7:end-5));
        
        dist = (CoM(i_img,2)-a)/cosd(CoM(i_img,3));
        img = imread([data_dirs{i} '/roi/roi_' num2str(i_img) '.png']);
        img_norm = imtranslate(img, [255 - CoM(i_img,1), 255 - CoM(i_img,2)]);
        img_norm = imrotate(img_norm, CoM(i_img,3));
        img_norm = imcrop(img_norm, [size(img_norm,1)/2-150 size(img_norm,2)/2-150 300 300]);
        
        img_norm2 = imtranslate(imori, [255 - CoM(i_img,1), 255 - CoM(i_img,2)]);
        img_norm2 = imrotate(img_norm2, CoM(i_img,3));
        img_norm2 = imcrop(img_norm2, [size(img_norm2,1)/2-150 size(img_norm2,2)/2-150 300 300]);
        
        [Y,X] = find(img_norm);
        length(j) = max(Y(find(X==150)))-min(Y(find(X==150)));
        ratio{i}(j+1) = abs(dist/length(j));
    end
    ratio{i}(1) =mean(ratio{i}(2:end));
    csvwrite([data_dirs{i} '/ratio.csv'], ratio{i}');
end

for i = 1:20
    mratio(i)= ratio{1,i}(1);
end
mean(mratio)
std(mratio)

fileID = fopen(['./Centroid_Position.csv'],'w');
new_ratio = zeros(30,20);
for i = 1:20
    fprintf(fileID,'%s \t', data_dirs{i}(12:end));
    for j = 1:size(ratio{i},2)
        new_ratio(j,i)=ratio{i}(j);
    end
end
fprintf(fileID,'\n');
fclose(fileID);
dlmwrite(['./Centroid_Position.csv'],new_ratio,'delimiter','\t','-append');