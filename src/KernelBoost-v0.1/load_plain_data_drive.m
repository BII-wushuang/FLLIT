function [imgs] = load_plain_data_drive(params,op,ch)

% Load plain image data
%
%  authors: Carlos Becker, Roberto Rigamonti, CVLab EPFL
%  e-mail: name <dot> surname <at> epfl <dot> ch
%  web: http://cvlab.epfl.ch/
%  date: February 2014

% imgs_no = params.(ch).(op).imgs_no;

opts = init_DRIVE;


switch op
    
    case 'train'
        
        ImgDir = [opts.dataDir opts.trnImgDir];
        
        imgIds = opts.imgIds;
        
        
        
       
        
    case 'test'
        
        ImgDir = [opts.dataDir opts.testImgDir];
        
        imgIds = opts.testimgIds;
       
       
    otherwise
        
        
end

imgs_no = length(imgIds);

for i_img = 1 : imgs_no
    
    imgs_list{i_img} = [ImgDir imgIds{i_img}];
    
end

imgs = cell(imgs_no,3);

for i_img = 1:imgs_no
    if strfind(imgs_list{i_img},'.nrrd')
        imgs{i_img} = permute(nrrdLoad(imgs_list{i_img}),[2,1]);
    else
        I = imread(imgs_list{i_img});
        
        for i_c =1 : size(I,3)
            
            imgs{i_img,i_c} = padarray(im2double(I(:,:,i_c)),[params.border_size,params.border_size],'symmetric');
            
        end

    end
    
    
    %imgs{i_img} = padarray(imgs{i_img},[params.border_size,params.border_size],'symmetric');
end

end
