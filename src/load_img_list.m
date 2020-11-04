function img_list = load_img_list(data_dir)
	formats = {'*.tif', '*.tiff', '*.bmp', '*.png', '*.jpg', '*.jpeg'};
	for i =1:length(formats)
		img_list = dir([data_dir formats{i}]);
		if(~isempty(img_list))
			break;
		end
	end
end