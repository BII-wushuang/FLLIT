function [imgs] = load_plain_data(params,op,ch)

% Load plain image data
%
%  authors: Carlos Becker, Roberto Rigamonti, CVLab EPFL
%  e-mail: name <dot> surname <at> epfl <dot> ch
%  web: http://cvlab.epfl.ch/
%  date: February 2014

imgs_no = params.(ch).(op).imgs_no;
imgs_list = params.(ch).(op).list;
imgs = cell(imgs_no,1);

for i_img = 1:imgs_no
    if strfind(imgs_list{i_img},'.nrrd')
        imgs{i_img} = permute(nrrdLoad(imgs_list{i_img}),[2,1]);
    else
        imgs{i_img} = im2double(imread(imgs_list{i_img}));
        if (size(imgs{i_img},3)>1)
            imgs{i_img} = rgb2gray(imgs{i_img});
        end
    end
    imgs{i_img} = padarray(imgs{i_img},[params.border_size,params.border_size],'symmetric');
end

end
