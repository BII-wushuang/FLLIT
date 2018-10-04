function [imgs] = load_binary_channel_L(params,op,ch)

% Load an image set as binary images
%
%  authors: Carlos Becker, Roberto Rigamonti, CVLab EPFL
%  e-mail: name <dot> surname <at> epfl <dot> ch
%  web: http://cvlab.epfl.ch/
%  date: February 2014


% switch params.dataset_name
%     
%     case 'neuron'
%                        
%         opts = init_Neuron_v4(params.opts);
% 


% switch ch
%     
%     case 'gts'
%         
%         imgIds = params.opts.imgIds;
%         
%         
%         
%         
% end

opts = params.opts;

    
    switch op
        
        case 'train'
            
            imgIds = params.opts.imgIds;
            
            gt_imgIds = params.opts.gt_imgIds;
            
            GtDir = opts.trnGtDir;
            
            imgDir = opts.trnImgDir;
            
        case 'test'
            
            imgIds = params.opts.testimgIds;
            
            gt_imgIds = params.opts.testgt_imgIds;
            
            GtDir = opts.testGtDir;
            
            imgDir = opts.testImgDir;
        otherwise
            
    end
    
       
    
    imgs_no = length(imgIds);
    
    imgs = cell(imgs_no,1);
    
    % imgs{i_img} = imread(params.(ch).(op).list{i_img});
    
    for i_img = 1:imgs_no
        
        
        switch params.dataset_name
            
            case 'road'
                
                 I = imread([opts.dataDir  GtDir gt_imgIds{i_img}]);
                 
                 
                 switch ch
                     
                     case 'gts'
                         
                         I = I > 255 / 3 * 2;
                         
                     case 'masks'
                         
                         I = ones(size(I));
                    
                     otherwise
                         
                 end
                         
            
            
            
            
            case 'neuron'
                
                
                
                
                I = imread([opts.dataDir  GtDir gt_imgIds{i_img}]);
                
                switch ch
                    
                    case 'gts'
                                                
                        I = I > 255 / 3 * 2;
                        
                    case 'masks'
                        
                        I = double(I < 255 / 3) + double(I > 255 / 3 * 2);
                        
                        I = I > 0;
                        
                    otherwise
                        
                end
                
                
            case 'STARE'
                
                
                
                switch ch
                    
                    case 'gts'
                        
                        I = imread([opts.dataDir  GtDir gt_imgIds{i_img}]);
                                                
                        I = I > opts.pos_thres;
                        
                    case 'masks'
                        
                        I = imread([opts.dataDir  imgDir imgIds{i_img}]);
                        
                        I = mean(I,3) > 40;
                        
%                         I = I > 0;
                        
                    otherwise
                        
                end               
                
            otherwise
                
        end
        
        if (~islogical(I))
            I = im2bw(I);
        end
        imgs{i_img} = padarray(I,[params.border_size,params.border_size],'symmetric');
        
        
        
    end
    
    
    
    
    
    
end








