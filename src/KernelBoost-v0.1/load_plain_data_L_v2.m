function [imgs] = load_plain_data_L_v2(params,op,ch)

% Load plain image data
%
%  authors: Carlos Becker, Roberto Rigamonti, CVLab EPFL
%  e-mail: name <dot> surname <at> epfl <dot> ch
%  web: http://cvlab.epfl.ch/
%  date: February 2014

opts = params.opts;

    switch op
        
        case 'train'
            
            imgIds = params.opts.imgIds;
            
            gt_imgIds = params.opts.gt_imgIds;
            
            ImgDir = opts.trnImgDir;
            
            GtDir = opts.trnGtDir;
            
            
        case 'test'
            
            imgIds = params.opts.testimgIds;
            
            gt_imgIds = params.opts.testgt_imgIds;
            
            ImgDir = opts.testImgDir;
            
            GtDir = opts.testGtDir;
            
        otherwise
            
    end
    
    


imgs_no = length(imgIds);
imgs_list = params.(ch).(op).list;
imgs = cell(imgs_no,1);

opts = params.opts;

for i_img = 1:imgs_no
    
    switch params.dataset_name
        
        case 'neuron'
                
            switch ch
                
                case 'imgs'
                                        
                    I = import_img([opts.dataDir  ImgDir imgIds{i_img}],'neuron');
                    
                    I = im2double(I(:,:,2));                    
                   
                otherwise
                    
            end
            
            
        case 'STARE'
            
        switch ch
            
            case 'imgs'
                
                I = imread([opts.dataDir  ImgDir imgIds{i_img}]);
                
                I = im2double(I(:,:,2)); 
                
                
            otherwise
                
        end
        
        case 'road'
            
            switch ch
                
                case 'imgs'
                    
                    I = imread([opts.dataDir  ImgDir imgIds{i_img}]);
                    
                    I = im2double(rgb2gray(I));
                    
                    
                otherwise
                    
            end
            
            
            
        otherwise
            
    end
    
    
%     if strfind(imgs_list{i_img},'.nrrd')
%         imgs{i_img} = permute(nrrdLoad(imgs_list{i_img}),[2,1]);
%     else
%         imgs{i_img} = im2double(imread(imgs_list{i_img}));
%         if (size(imgs{i_img},3)>1)
%             imgs{i_img} = rgb2gray(imgs{i_img});
%         end
%     end
    imgs{i_img} = padarray(I,[params.border_size,params.border_size],'symmetric');
end

end
