function [imgs] = load_binary_channel(params,op,ch)

% Load an image set as binary images
%
%  authors: Carlos Becker, Roberto Rigamonti, CVLab EPFL
%  e-mail: name <dot> surname <at> epfl <dot> ch
%  web: http://cvlab.epfl.ch/
%  date: February 2014

imgs_no = params.(ch).(op).imgs_no;
imgs = cell(imgs_no,1);

for i_img = 1:imgs_no
    imgs{i_img} = imread(params.(ch).(op).list{i_img});
    if (~islogical(imgs{i_img}))
        imgs{i_img} = im2bw(imgs{i_img});
    end
    imgs{i_img} = padarray(imgs{i_img},[params.border_size,params.border_size],'symmetric');
    
end

end
