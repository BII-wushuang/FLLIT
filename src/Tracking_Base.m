% A very simple tracker based on the hungarian linker method
% Outputs tip positions and identities of individual legs
function Tracking_Base (data_dir)
nLegs = 6; %preassume there are 6 legs

%% Section 1: locate image folder and create output folder
if (nargin < 1)
    data_dir = uigetdir('./Data');
end

pos_bs = strfind(data_dir,'Data');
sub_dir = data_dir(pos_bs(end)+length('Data'):length(data_dir));
seg_dir = ['./Results/SegmentedImages' sub_dir '/'];
data_dir = ['./Data' sub_dir '/'];
img_list = load_img_list(data_dir);

output_dir = ['./Results/Tracking' sub_dir '/'];
if(~exist(output_dir))
    mkdir(output_dir);
end
%% End of section 1
%% Section 2:
se = strel('disk',4);

% maximal movement of legs in between frames
max_distance = 20;

trajectory = [];
norm_trajectory = [];
CoM = [];

end_frame = length(dir([seg_dir '*.png']))/2;
for i = 1: end_frame
    clear connComp;
    % Load necessary images
    Im = (imread([seg_dir 'img_' num2str(i) '.png']));
    Imroi = imread([seg_dir 'roi_' num2str(i) '.png']);

    % Fix image/segmentation irregularities
    Imwork = double(Im(:,:,1) == 255) .* Imroi;
    connComp = bwconncomp(Imwork);

    if (connComp.NumObjects == 6 && min(cellfun(@length,connComp.PixelIdxList)) > 35)
        Imroi = imerode(Imroi,se);
        [locY,locX] = find(Imroi > 0);

        [centre_of_mass,theta] = findCoM(locX,locY);
        locX = locX - centre_of_mass(1);
        locY = locY - centre_of_mass(2);
        points = [cosd(theta) sind(theta); -sind(theta) cosd(theta)]* [locX'; locY'];

        %Attempt to ensure that the head is in the positive
        %y-direction by comparing which side has more number of
        %pixel points
        if (length(find(points(2,:) > 0))< length(find(points(2,:) < 0)))
            theta = mod(theta + 180,360);
        end
        
        for cc = 1:nLegs
            legpixels = connComp.PixelIdxList{cc};

            imleg = zeros(size(Imwork));
            imleg(legpixels) = 1;

            [skr,~] = skeleton(imleg);
            skr_start = 3;
            [~,exy,~] = anaskel(bwmorph(skr > skr_start,'skel',Inf));
%             while (length(exy)>2)
%                 skr_start = skr_start+1;
%                 [~,exy,~] = anaskel(bwmorph(skr > skr_start,'skel',Inf));
%             end

            imgX = exy(1,:) - centre_of_mass(1);
            imgY = exy(2,:) - centre_of_mass(2);
            normpoints = [cosd(theta) sind(theta); -sind(theta) cosd(theta)]* [imgX; imgY];

            [~,tipidx]=max(min(pdist2(points',normpoints','euclidean')));

            for j = 1:2
                raw_normtips(cc,j) = normpoints(j,tipidx);
                raw_tips(cc,j) = exy(j,tipidx);
            end
        end

        leftlegs_idx = find(raw_normtips(:,1)<0);
        if(length(leftlegs_idx)~=3)
            continue;
        end
        rightlegs_idx = find(raw_normtips(:,1)>=0);

        [~,left_sort] = sort(raw_normtips(leftlegs_idx,2));

        leg_counter = 0;
        for kk = 1 :3
            trajectory(i, kk, :) = raw_tips(leftlegs_idx(left_sort(kk-leg_counter)),:);
            norm_trajectory(i, kk, :) = raw_normtips(leftlegs_idx(left_sort(kk-leg_counter)),:);
        end

        leg_counter = 0;
        [~,right_sort] = sort(raw_normtips(rightlegs_idx,2));
        for kk = 1 :3
            trajectory(i, kk+3, :) = raw_tips(rightlegs_idx(right_sort(kk-leg_counter)),:);
            norm_trajectory(i, kk+3, :) = raw_normtips(rightlegs_idx(right_sort(kk-leg_counter)),:);
        end

        CoM(i, :) = [centre_of_mass, theta];

        start_frame = i;
        fprintf('Tracking automatically initiated on frame %d for dataset %s\n', i, sub_dir);
        for k = i-1 : -1 : 1
            Im = (imread([seg_dir 'img_' num2str(k) '.png']));
            Imroi = imread([seg_dir 'roi_' num2str(k) '.png']);

            % Fix image/segmentation irregularities
            Imwork = double(Im(:,:,1) == 255) .* Imroi;
            Imwork = Imwork - (bwdist(imerode(Imroi,se))<10);
            [locY,locX] = find(imerode(Imroi,se) > 0);
            
            %Centre_of_Mass
            [centre_of_mass,theta] = findCoM(locX,locY);

            if (abs(theta - CoM(k+1,3)) > 90 && abs(theta - CoM(k+1,3)) < 270)
                theta = mod(theta + 180,360);
            end

            CoM(k, :) = [centre_of_mass, theta];

            [raw_tips, raw_normtips] = findtips(Imwork, Imroi, locX, locY, centre_of_mass, theta);
            
            x_t1 = raw_normtips;
            x_t0 = reshape(norm_trajectory(k+1,:),[nLegs 2]);
            target_indices = hungarianlinker(x_t0, x_t1, max_distance);

            for kk = 1:length(target_indices)
                if (target_indices(kk) > 0)
                    norm_trajectory(k, kk, :) = raw_normtips(target_indices(kk),:);
                    trajectory(k, kk, :) = raw_tips(target_indices(kk),:);
                else
                    norm_trajectory(k, kk, :) = 0;
                    trajectory(k, kk, :) = 0;
                end
            end

            leftoveridx = find(hungarianlinker(raw_normtips, x_t0, max_distance)==-1);
            x_t1 = raw_normtips(leftoveridx,:);
            x_t2 = raw_tips(leftoveridx,:);

            last_seen_tips = reshape(norm_trajectory(k,:),[nLegs 2]);
            unassigned_idx = find(last_seen_tips(:,1) == 0);
            for kk = 1:length(unassigned_idx)
                idx = find(norm_trajectory(1:i,unassigned_idx(kk),1),1,'last');
                last_seen_tips(unassigned_idx(kk),:) = reshape(norm_trajectory(idx,unassigned_idx(kk),:), [1 2]);
            end

            if (~isempty(unassigned_idx) && ~isempty(x_t1))
                    C = hungarianlinker(last_seen_tips(unassigned_idx,:), x_t1, 1.25*max_distance);
                    for l = 1:length(C)
                        if (C(l) ~= -1)
                            norm_trajectory(k,unassigned_idx(l), :) = x_t1(C(l), :);
                            trajectory(k,unassigned_idx(l), :) = x_t2(C(l), :);
                        end
                    end
            end
        end
        break
    else
        continue;
    end
end

skip = 1;

for i = start_frame+1 : skip : end_frame
    if (mod(i,50) == 1)
        fprintf('Tracking progress: frame %d for dataset %s\n', i, sub_dir);
    end
    % Load necessary images
    Im = (imread([seg_dir 'img_' num2str(i) '.png']));
    Imroi = imread([seg_dir 'roi_' num2str(i) '.png']);

    % Fix image/segmentation irregularities
    Imwork = double(Im(:,:,1) == 255) .* Imroi;
    Imwork = Imwork - (bwdist(imerode(Imroi,se))<10);
    [locY,locX] = find(imerode(Imroi,se) > 0);

    %Centre_of_Mass
    [centre_of_mass,theta] = findCoM(locX,locY);

    if (abs(theta - CoM(i-1,3)) > 90 && abs(theta - CoM(i-1,3)) < 270)
        theta = mod(theta + 180,360);
    end

    CoM(i, :) = [centre_of_mass, theta];

    [raw_tips, raw_normtips] = findtips(Imwork, Imroi, locX, locY, centre_of_mass, theta);
    
    x_t1 = raw_normtips;
    x_t0 = reshape(norm_trajectory(i-1,:),[nLegs 2]);
    target_indices = hungarianlinker(x_t0, x_t1, max_distance);

    for kk = 1:length(target_indices)
        if (target_indices(kk) > 0)
        norm_trajectory(i, kk, :) = raw_normtips(target_indices(kk),:);
        trajectory(i, kk, :) = raw_tips(target_indices(kk),:);
        end
    end

    leftoveridx = find(hungarianlinker(raw_normtips, x_t0, max_distance)==-1);
    x_t1 = raw_normtips(leftoveridx,:);
    x_t2 = raw_tips(leftoveridx,:);

    last_seen_tips = reshape(norm_trajectory(i,:),[6 2]);
    unassigned_idx = find(last_seen_tips(:,1) == 0);
    for kk = 1:length(unassigned_idx)
        idx = find(norm_trajectory(:,unassigned_idx(kk),1),1,'last');
        last_seen_tips(unassigned_idx(kk),:) = reshape(norm_trajectory(idx,unassigned_idx(kk),:), [1 2]);
    end

    if (~isempty(unassigned_idx) && ~isempty(x_t1))
        C = hungarianlinker(last_seen_tips(unassigned_idx,:), x_t1, max_distance);
        for l = 1:length(C)
            if (C(l) ~= -1)
                norm_trajectory(i,unassigned_idx(l), :) = x_t1(C(l), :);
                trajectory(i,unassigned_idx(l), :) = x_t2(C(l), :);
            end
        end
    end
end

save([output_dir 'CoM.mat'],'CoM');
save([output_dir 'trajectory.mat'],'trajectory');
save([output_dir 'norm_trajectory.mat'],'norm_trajectory');
csvwrite([output_dir 'CoM.csv'],CoM);
csvwrite([output_dir 'trajectory.csv'],trajectory);
csvwrite([output_dir 'norm_trajectory.csv'],norm_trajectory);
end