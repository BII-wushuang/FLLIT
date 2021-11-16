% Processes the raw data to extract relevant parameters
function DataProcessing_New(fps,data_dir,scale, bodylength)

if (nargin < 1)
    fps = 2000;
    scale=11; % dimension of arena is 11mm
    data_dir = uigetdir('./Results/Tracking/');
    bodylength = 2.88;
end

bodylength = bodylength * 512 / scale;

try
    addpath('./Export-Fig');
catch
end

pixel_mm = 11/512; % 1 pixel corresponds to how many mm
frame_ms = 1000/fps; % 1 frame corresponds to how many ms

pos_bs = strfind(data_dir,'Results');
sub_dir = data_dir(pos_bs(end)+length('Results/Tracking/'):length(data_dir));
output_dir = ['./Results/Tracking/' sub_dir '/'];
seg_dir = ['./Results/SegmentedImages/' sub_dir '/'];

load([output_dir 'trajectory.mat']);
load([output_dir 'norm_trajectory.mat']);
load([output_dir 'CoM.mat']);

%trajectory(:,:,2) = -trajectory(:,:,2);
norm_trajectory(:,:,2) = -norm_trajectory(:,:,2);

nlegs = size(trajectory,2);
if (nlegs == 6)
    legs_id = {'L1', 'L2', 'L3', 'R1', 'R2', 'R3'};
    colors = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 0 1 1; 1 0 1];
else
    legs_id = {'L1', 'L2', 'L3', 'L4', 'R1', 'R2', 'R3', 'R4'};
    colors = [1 0 0; 0 1 0; 0 0 1; 1 0 0.5; 1 1 0; 1 0 1; 0 1 1; 1 0.5 0];
end

start = zeros(nlegs,1);
missing = zeros(nlegs,1);
for j = 1:nlegs
    if (isempty(find(trajectory(:,j,1),1)))
        missing(j) = 1;
        continue;
    end
    start(j) = find(trajectory(:,j,1),1);
end
startframe = max(start);
endframe = length(CoM);
nframes = endframe - startframe + 1;
dist = zeros(endframe,1);
x = zeros(endframe, nlegs);
y = zeros(endframe, nlegs);
nx = zeros(endframe, nlegs);
ny = zeros(endframe, nlegs);

% body length
try
    img = imread([seg_dir 'roi_' num2str(1) '.png']);
    
    body_length = zeros(nframes,1);
    for i = startframe : endframe
        img = imread([seg_dir 'roi_' num2str(i) '.png']);
        img_norm = imtranslate(img, [255 - CoM(i,1), 255 - CoM(i,2)]);
        img_norm = imrotate(img_norm, CoM(i,3));
        img_norm = imcrop(img_norm, [size(img_norm,1)/2-150 size(img_norm,2)/2-150 300 300]);
        [Y,X] = find(img_norm);
        body_length(i-startframe+1) = max(Y(find(X==150)))-min(Y(find(X==150)));
    end
    fileID = fopen([output_dir 'bodylength.xlsx'],'w');
    fprintf(fileID, '%s\t\n', 'Mean body length (mm)');
    fprintf(fileID, '%d\n\n', mean(body_length)*scale/512);
    fprintf(fileID, '%s\t%s\t\n', 'Frame #', 'Length (mm)');
    fclose(fileID);
    dlmwrite([output_dir 'bodylength.xlsx'],[(startframe:endframe)', body_length*scale/512],'delimiter','\t','-append');
catch
end

for j =1:nlegs
    if(missing(j))
        continue;
    end
    nonzero = find(trajectory(:,j,1));
    for i = startframe : endframe
        if (trajectory(i,j,1) ~= 0)
            x(i,j) = trajectory(i,j,1);
            y(i,j) = trajectory(i,j,2);
            nx(i,j) = norm_trajectory(i,j,1);
            ny(i,j) = norm_trajectory(i,j,2);
        else
            previdx = nonzero(find(nonzero<i, 1));
            nextidx = nonzero(find(nonzero>i, 1));
            if (~isempty(nextidx))
                x(i,j) = trajectory(previdx,j,1) + (i-previdx)/(nextidx-previdx) * (trajectory(nextidx,j,1) - trajectory(previdx,j,1));
                y(i,j) = trajectory(previdx,j,2) + (i-previdx)/(nextidx-previdx) * (trajectory(nextidx,j,2) - trajectory(previdx,j,2));
                nx(i,j) = norm_trajectory(previdx,j,1) + (i-previdx)/(nextidx-previdx) * (norm_trajectory(nextidx,j,1) - norm_trajectory(previdx,j,1));
                ny(i,j) = norm_trajectory(previdx,j,2) + (i-previdx)/(nextidx-previdx) * (norm_trajectory(nextidx,j,2) - norm_trajectory(previdx,j,2));
            else
                x(i,j) = x(i-1,j);
                y(i,j) = y(i-1,j);
                nx(i,j) = nx(i-1,j);
                ny(i,j) = ny(i-1,j);
            end
       end
    end
end

for i = startframe : endframe
    if (i == startframe)
        dist(i) = 0;
    else
        dist(i) = dist(i-1) + (CoM(i,1) - CoM(i-1,1))*sind(CoM(i-1,3)) + (CoM(i,2) - CoM(i-1,2))*-cosd(CoM(i-1,3));
    end
end

%% Locate turning points
%Finding the turning points: first use the RDP algorithm to reduce number
%of points to describe the 2D trajectory, if the turning angle at a reduced point 
%is larger than a threshold (set to 50 deg), this point
%will be identified as a turning point.
reducedpoints = DouglasPeucker([CoM(startframe:endframe,1), CoM(startframe:endframe,2)],10,0);
turn_angle = zeros(length(reducedpoints)-2,1);
for i = 1 : length(reducedpoints) - 2
    vec_1 = [reducedpoints(i+2,1)-reducedpoints(i+1,1), reducedpoints(i+2,2)-reducedpoints(i+1,2)];
    vec_2 = [reducedpoints(i+1,1)-reducedpoints(i,1), reducedpoints(i+1,2)-reducedpoints(i,2)];
    turn_angle(i) = acosd(dot(vec_1,vec_2) / (norm(vec_1)*norm(vec_2)));    
end
turn_thres = 50;
turn_points_idx = find(turn_angle>turn_thres);
turn_points = zeros(length(turn_points_idx),2);

for i = 1 : length(turn_points_idx)
    turn_points(i,:) = [reducedpoints(turn_points_idx(i)+1,1), reducedpoints(turn_points_idx(i)+1,2)];
    turn_labels{i} = ['Turn' num2str(i) ' = (' num2str(round(turn_points(i,1))) ', ' num2str(round(turn_points(i,2))) ')'];
end

try
    clf(1);
catch
end
figure(1);
scatter(CoM(startframe : endframe,1),CoM(startframe : endframe,2),'*');
title('Trajectory');
axis([0 512 0 512]);
axis square
set(gca,'Ydir','reverse')
% if (~isempty(turn_points_idx))
%     labelpoints(turn_points(:,1)', turn_points(:,2)', turn_labels, 'Center');
% end
export_fig([output_dir 'BodyTrajectory.pdf'], '-pdf','-transparent');

%% Obtaining the CoM speed
i = startframe : endframe;
Distance = dist(startframe : endframe);
skip = min(nframes/2,floor(50/frame_ms));
B(1) = 0;
B(startframe+1:endframe) = diff(dist(startframe:endframe));
f = fit((startframe + skip / 2: skip: endframe - skip/2)', B(startframe + skip / 2: skip: endframe - skip/2)', 'smoothingspline', 'SmoothingParam', 0.8);
%Speedhist = hist(11000/512*f(i),(-20:0.5:50));
Speed = [[1:nframes]'*frame_ms, f(i)*pixel_mm*1000/frame_ms];

fileID = fopen([output_dir 'BodyVelocity.xlsx'],'w');
fprintf(fileID,'%s \t %s', 'Time(ms)', 'Velocity (mm/s)');
fprintf(fileID,'\n');
fclose(fileID);
dlmwrite([output_dir 'BodyVelocity.xlsx'],Speed,'delimiter','\t','-append');

clf(1);
figure(1);
hold on
Speedplot = plot(Speed(:,1),Speed(:,2));
plot(zeros(nframes,1),'--','color',[0 0 0]);
hold off
title('Body Velocity');
axis([1 nframes*frame_ms -30 60]);
Speedplot.Parent.XLabel.String = 'Time(ms)';
Speedplot.Parent.YLabel.String = 'Velocity (mm/s)';
export_fig([output_dir 'BodyVelocity.pdf'], '-pdf','-transparent');

% Speedhistplot = plot(Speedhist);
% title('Histogram of Drosophila Velocity');
% Speedhistplot.Parent.XLabel.String = 'Velocity';
% Speedhistplot.Parent.YLabel.String = '# Occurence';
% export_fig([output_dir 'HistogramVelocity.pdf'], '-pdf','-transparent');

%% Leg Speed and Gait
legspd = zeros(nframes,nlegs);
nlegspd = zeros(nframes,nlegs);
legsy = zeros(nframes,nlegs);
gait = zeros(nframes,nlegs);
ngait = zeros(nframes,nlegs);
foot_dragging = zeros(nframes,nlegs);
i = startframe : endframe;

clf(1);
figure(1);
for j = 1:nlegs
    % Speed in arena-centered frame of reference
    displacement(startframe: endframe)=sqrt(x(startframe: endframe,j).^2+y(startframe: endframe,j).^2);
    for k = startframe:endframe
        if (k==startframe)
            spd(k,j) = 0;
        else
            spd(k,j) = sqrt((x(k,j)-x(k-1,j))^2+(y(k,j)-y(k-1,j))^2);
        end
    end
    f = fit((startframe: 10: endframe)', spd(startframe: 10: endframe,j), 'smoothingspline', 'SmoothingParam', 1);
    legspd(:, j) = abs(f(i))*pixel_mm*1000/frame_ms;
    
    % Vertical leg trajectories in body-centred frame of reference
    f = fit((startframe: 10: endframe)', ny(startframe: 10: endframe,j), 'smoothingspline', 'SmoothingParam', 0.95);
    legsy(:, j) = f(i);
    
    hold on
    subplot(2, round(nlegs/2), j);
    plot(legsy(:,j));
    title(legs_id{j});
    axis([1 nframes -150 150]);
    hold off
    
    % Gait coordination
    skip = floor(10 / frame_ms);
    g = fit((startframe: skip: endframe)', x(startframe: skip: endframe,j), 'smoothingspline', 'SmoothingParam', 1);
    h = fit((startframe: skip: endframe)', y(startframe: skip: endframe,j), 'smoothingspline', 'SmoothingParam', 1);
    gait_idx = find(abs(gradient(g(i)))>0.3*frame_ms | abs(gradient(h(i)))>0.3*frame_ms);
    % gait_idx = find(legspd(:,j)>mean(legspd(:,j))/3);
    gait(gait_idx,j) = 1;
    
    ng = fit((startframe: skip: endframe)', nx(startframe: skip: endframe,j), 'smoothingspline', 'SmoothingParam', 1);
    nh = fit((startframe: skip: endframe)', ny(startframe: skip: endframe,j), 'smoothingspline', 'SmoothingParam', 1);
    ngait_idx = find(abs(gradient(ng(i)))>0.2*frame_ms | abs(gradient(nh(i)))>0.2*frame_ms);
    ngait(ngait_idx,j) = 1;
    
    thres = 15/frame_ms;
    cc = bwconncomp(gait(:,j));
    for k = 1 : cc.NumObjects
        if cellfun(@length,cc.PixelIdxList(k))<thres
            gait(cc.PixelIdxList{k},j) = 0;
        end
    end
    
    cc = bwconncomp(~gait(:,j));
    for k = 1 : cc.NumObjects
        if cellfun(@length,cc.PixelIdxList(k))<thres
            gait(cc.PixelIdxList{k},j) = 1;
        end
    end
    
    cc = bwconncomp(ngait(:,j));
    for k = 1 : cc.NumObjects
        if cellfun(@length,cc.PixelIdxList(k))<thres
            ngait(cc.PixelIdxList{k},j) = 0;
        end
    end
    
    cc = bwconncomp(~ngait(:,j));
    for k = 1 : cc.NumObjects
        if cellfun(@length,cc.PixelIdxList(k))<thres
            ngait(cc.PixelIdxList{k},j) = 1;
        end
    end
    
    % Foot dragging is identified as an event where the leg is moving with
    % respect to the arena yet relatively stationary with respect to the
    % body. The episode must last a minimum of 20ms.
    foot_dragging(:,j) = ~ngait(:,j).*gait(:,j);
    cc = bwconncomp(foot_dragging(:,j));
    for k = 1 : cc.NumObjects
        if cellfun(@length,cc.PixelIdxList(k))<20/frame_ms
            foot_dragging(cc.PixelIdxList{k},j) = 0;
        end
    end
    
    clear connComp;
end
mtit('Vertical leg trajectories in body-centred frame of reference','xoff',0,'yoff',.04);
export_fig([output_dir 'LegVerticalTrajectories.pdf'], '-pdf','-transparent');

% legspd = legspd.*gait;

fileID = fopen([output_dir 'LegSpeed.xlsx'],'w');
fprintf(fileID,'%s \t %s \t', 'Time(ms)', legs_id{:});
fprintf(fileID,'\n');
fclose(fileID);
dlmwrite([output_dir 'LegSpeed.xlsx'],[[1:nframes]'*frame_ms, legspd],'delimiter','\t','-append');

% fileID = fopen([output_dir 'FootDragging.xlsx'],'w');
% fprintf(fileID,'%s \t %s \t %s \n', 'Leg', 'StartFrame', 'EndFrame');
% for j = 1 : nlegs
%     cc = bwconncomp(foot_dragging(:,j));
%     if(cc.NumObjects>0)
%         fprintf(fileID,'%s', legs_id{j});
%         for k = 1 : cc.NumObjects
%             fprintf(fileID,'\t %d \t %d \n', min(cc.PixelIdxList{k}), max(cc.PixelIdxList{k}));
%         end
%     end
% end
% fclose(fileID);

clf(1);
figure(1);
colormap('jet');
clims = [0 200];
Gaitplot = imagesc([1 nframes*frame_ms], [1 nlegs], abs(legspd)', clims);
colorbar;
title('Gait and Speed');
Gaitplot.Parent.YTickLabel = legs_id;
Gaitplot.Parent.XLabel.String = 'Time(ms)';
export_fig([output_dir 'Gait.pdf'], '-pdf','-transparent');

%% Gait index
if (nlegs ==6)
gaitindex = zeros(nframes,1);
for i = 1:nframes
    if  isequal(gait(i,:),[1 0 1 0 1 0]) ||  isequal(gait(i,:),[0 1 0 1 0 1]) 
        gaitindex(i) = 1;
        continue;
    elseif  isequal(gait(i,:),[1 0 0 0 1 0]) ||  isequal(gait(i,:),[0 1 0 1 0 0 0]) ||  isequal(gait(i,:),[1 0 0 0 0 1]) ||  isequal(gait(i,:),[0 0 1 1 0 0]) ||  isequal(gait(i,:),[0 1 0 0 0 1]) ||  isequal(gait(i,:),[0 0 1 0 1 0])
        gaitindex(i) = -1;
    end
end
% Gait index averaged over a moving window of 80ms
gaitindex = movmean(gaitindex,80/frame_ms);

fileID = fopen([output_dir 'GaitIndex.xlsx'],'w');
fprintf(fileID,'%s \t %s', 'Time(ms)', 'Gait Index');
fprintf(fileID,'\n');
fclose(fileID);
dlmwrite([output_dir 'GaitIndex.xlsx'],[[1:nframes]'*frame_ms, gaitindex],'delimiter','\t','-append');

clf(1);
figure(1);
GaitIndexplot = plot([1:nframes]'*frame_ms, gaitindex);
title('Gait Index');
axis([0 nframes*frame_ms -1 1]);
GaitIndexplot.Parent.YLabel.String = 'Gait Index';
GaitIndexplot.Parent.XLabel.String = 'Time(ms)';
export_fig([output_dir 'GaitIndex.pdf'], '-pdf','-transparent');
end

%% Stride Parameters
fileID = fopen([output_dir 'StrideParameters.xlsx'],'w');
fclose(fileID);

landing_std = zeros(nlegs,2);
takeoff_std = zeros(nlegs,2);
landing_mean = zeros(nlegs,2);
takeoff_mean = zeros(nlegs,2);
for j = 1:nlegs
    if(missing(j))
        continue;
    end
    cc = bwconncomp(gait(:,j));
    nstrides(j) = cc.NumObjects;
    
    if(nstrides(j)==0)
         fileID = fopen([output_dir 'StrideParameters.xlsx'],'a');
        fprintf(fileID,'%s',legs_id{j});
        fprintf(fileID,'\n');
        fprintf(fileID,'%s \t', 'Stride #', 'Duration (ms)', 'Period (ms)', 'Displacement (mm)', 'Path Covered (mm)', 'Take-off time (ms)', 'Landing time (ms)', 'AEP x (mm)', 'AEP y (mm)', 'PEP x (mm)', 'PEP y (mm)', 'Amplitude (mm)', 'Stance linearity (mm)', 'Stretch (mm)');
        fprintf(fileID,'\n');
        fclose(fileID);
        continue;
    end

    for k = 1 : nstrides(j)
        stride_start(j,k) = min(cc.PixelIdxList{k});
        stride_duration(j,k) = max(cc.PixelIdxList{k}) - min(cc.PixelIdxList{k});
        stride_end(j,k) = max(cc.PixelIdxList{k});
        stride_dist(j,k) = sqrt((y(startframe-1 + stride_start(j,k)+stride_duration(j,k),j) - y(startframe-1 + stride_start(j,k),j))^2+(x(startframe-1 + stride_start(j,k)+stride_duration(j,k),j) - x(startframe-1 + stride_start(j,k),j))^2);
        path_length(j,k) = sum(spd(startframe-1 + stride_start(j,k):startframe-1 + stride_start(j,k)+stride_duration(j,k),j));
        landing(j,k,:) = [nx(startframe-1 + stride_start(j,k)+stride_duration(j,k),j), ny(startframe-1 + stride_start(j,k)+stride_duration(j,k),j)];
        takeoff(j,k,:) = [nx(startframe-1 + stride_start(j,k),j), ny(startframe-1 + stride_start(j,k),j)];
        Xinter= interp1(startframe-1 + stride_start(j,k):min(20,stride_duration(j,k)-1):startframe-1 + stride_start(j,k) + stride_duration(j,k),nx(startframe-1 + stride_start(j,k):min(20,stride_duration(j,k)-1):startframe-1 + stride_start(j,k)+stride_duration(j,k),j),startframe-1+stride_start(j,k):startframe-1+stride_start(j,k)+stride_duration(j,k),'spline');
        Yinter= interp1(startframe-1 + stride_start(j,k):min(20,stride_duration(j,k)-1):startframe-1 + stride_start(j,k) + stride_duration(j,k),ny(startframe-1 + stride_start(j,k):min(20,stride_duration(j,k)-1):startframe-1 + stride_start(j,k)+stride_duration(j,k),j),startframe-1+stride_start(j,k):startframe-1+stride_start(j,k)+stride_duration(j,k),'spline');
        stride_amplitude(j,k) = ny(startframe-1 + stride_start(j,k)+stride_duration(j,k),j) - ny(startframe-1 + stride_start(j,k),j);
        stride_regularity(j,k) = mean(sqrt(sum(abs([nx(startframe-1+stride_start(j,k):startframe-1+stride_start(j,k)+stride_duration(j,k),j), ny(startframe-1+stride_start(j,k):startframe-1+stride_start(j,k)+stride_duration(j,k),j)] - [Xinter' Yinter']).^2,2)));
        stretch(j,k) =  mean(sqrt(sum(nx(startframe-1+stride_start(j,k)+round(stride_duration(j,k)/2),j)^2+(ny(startframe-1+stride_start(j,k)+round(stride_duration(j,k)/2),j)-0.157*bodylength)^2)));
    end
    
    for k = 1 : nstrides(j)-1
        stride_period(j,k) = stride_start(j,k+1) - stride_start(j,k);
    end
    if (nstrides(j)>1)
        mean_stride_period(j) = sum(stride_period(j,:))/(nstrides(j)-1);
    else
         stride_period(j,1) = nan;
         mean_stride_period(j) = nan;
    end
    stride_period(j,nstrides(j)) = nan;
    landing_std(j,:) = squeeze(std(landing(j,1 : nstrides(j),:)))';
    takeoff_std(j,:) = squeeze(std(takeoff(j,1 : nstrides(j),:)))';
    landing_mean(j,:) = squeeze(sum(landing(j,:,:)))'/nstrides(j);
    takeoff_mean(j,:) = squeeze(sum(takeoff(j,:,:)))'/nstrides(j);
    l_std(j) = 0;
    t_std(j) = 0;
    for k = 1: nstrides(j)
        l_std(j) = l_std(j) + sqrt((landing(j,k,1)-landing_mean(j,1))^2 + (landing(j,k,2)-landing_mean(j,2))^2);
        t_std(j) = t_std(j) + sqrt((takeoff(j,k,1)-takeoff_mean(j,1))^2 + (takeoff(j,k,2)-takeoff_mean(j,2))^2);
    end
    l_std(j) = l_std(j) / nstrides(j);
    t_std(j) = t_std(j) / nstrides(j);
    
    total_path_length(j) = sum(path_length(j,:));
    mov_percentage(j) = sum(stride_duration(j,:))/nframes;

    fileID = fopen([output_dir 'StrideParameters.xlsx'],'a');
    fprintf(fileID,'%s',legs_id{j});
    fprintf(fileID,'\n');
    fprintf(fileID,'%s \t', 'Stride #', 'Duration (ms)', 'Period (ms)', 'Displacement (mm)', 'Path Covered (mm)', 'Take-off time (ms)', 'Landing time (ms)', 'AEP x (mm)', 'AEP y (mm)', 'PEP x (mm)', 'PEP y (mm)', 'Amplitude (mm)', 'Stance linearity (mm)', 'Stretch (mm)');
    fprintf(fileID,'\n');
    fclose(fileID);
    dlmwrite([output_dir 'StrideParameters.xlsx'], [[1:nstrides(j)]' stride_duration(j,1:nstrides(j))'*frame_ms stride_period(j,1:nstrides(j))'*frame_ms stride_dist(j,1:nstrides(j))'*scale/512 path_length(j,1:nstrides(j))'*scale/512 stride_start(j,1:nstrides(j))'*frame_ms stride_end(j,1:nstrides(j))'*frame_ms landing(j,1:nstrides(j),1)'*scale/512 landing(j,1:nstrides(j),2)'*scale/512 takeoff(j,1:nstrides(j),1)'*scale/512 takeoff(j,1:nstrides(j),2)'*scale/512 stride_amplitude(j,1:nstrides(j))'*scale/512 stride_regularity(j,1:nstrides(j))'*scale/512 stretch(j,1:nstrides(j))'*scale/512],'delimiter','\t','-append');
    fileID = fopen([output_dir 'StrideParameters.xlsx'],'a');
    fprintf(fileID,'\n');
    fclose(fileID);
    dlmwrite([output_dir 'StrideParameters.xlsx'], ['______________'],'delimiter','\t','-append');
end

%% Path area covered
% Returns the area of the minimum convex hull required to cover the entire
% path covered by each leg in the body-centered frame of reference
area = zeros(nlegs,1);
area_over_path = zeros(nlegs,1);
pathpca = zeros(nlegs,2,2);
centre = zeros(nlegs, 2);

fileID = fopen([output_dir 'LegParameters.xlsx'],'w');
% fprintf(fileID,'%s \t', 'Leg', 'Movement %', 'Mean Stride Period (ms)', 'Total Path Covered (mm)', 'Mean Landing x (mm)', 'Mean Landing y (mm)', 'Mean Take-off x (mm)', 'Mean Take-off y (mm)', 'Std Landing x (mm)', 'Std Landing y (mm)', 'Std Take-off x (mm)', 'Std Take-off y (mm)', 'Domain Area (mm^2)', 'Domain Area / Path (mm)', 'Domain Length (mm)', 'Domain Width (mm)', 'Footprint/PCA deviation (deg)');
fprintf(fileID,'%s \t', 'Leg', 'Movement %', 'Mean Stride Period (ms)', 'Total Path Covered (mm)', 'Mean AEP x (mm)', 'Mean AEP y (mm)', 'Mean PEP x (mm)', 'Mean PEP y (mm)', 'Std AEP (mm)', 'Std PEP (mm)', 'Domain Area (mm^2)', 'Domain Area / Path (mm)', 'Domain Length (mm)', 'Domain Width (mm)', 'Footprint/PCA deviation (deg)');
fprintf(fileID,'\n');
clf(1);
figure(1);
for j = 1:nlegs
    if(missing(j) || nstrides(j)==0)
        continue;
    end
    [convex{j},area(j)] = convhull(nx(startframe : endframe,j),ny(startframe : endframe,j));
    area_over_path(j) = area(j)/total_path_length(j);
    centre(j,:) = mean([nx(startframe : endframe,j),ny(startframe : endframe,j)]);
    [domain_length(j), domain_width(j)] = projectpca([nx(startframe : endframe,j),ny(startframe : endframe,j)]);
    
    bounding = boundary(nx(startframe : endframe,j),ny(startframe : endframe,j));
    bounding = mod(bounding,nframes)+((mod(bounding,nframes)==0)*nframes)+startframe-1;
    hold on
    legplot(j) = scatter(nx(startframe : endframe,j),ny(startframe : endframe,j),3,[colors(j,1) colors(j,2) colors(j,3)],'filled');
    landingplot(j) = scatter(landing_mean(j,1),landing_mean(j,2),50,[0 0 0],'filled');
    takeoffplot(j) = scatter(takeoff_mean(j,1),takeoff_mean(j,2),50,[0 0 0],'d','filled');
    landingplot(j) = scatter(landing_mean(j,1),landing_mean(j,2),30,[colors(j,1) colors(j,2) colors(j,3)],'filled');
    takeoffplot(j) = scatter(takeoff_mean(j,1),takeoff_mean(j,2),30,[colors(j,1) colors(j,2) colors(j,3)],'d','filled');
    % plot(nx(bounding,j),ny(bounding,j));
    k = -100:100;
   
    stride_vec = [landing_mean(j,1)-takeoff_mean(j,1) landing_mean(j,2)-takeoff_mean(j,2)];
    stride_vec = stride_vec/norm(stride_vec);
    footprint_dev(j) = acos(stride_vec(1)*pathpca(j,1,1) + stride_vec(2)*pathpca(j,1,2))*180/pi;
    hold off
    fprintf(fileID,'%s',legs_id{j});
    fprintf(fileID,'\t %0.2f', [100*mov_percentage(j) mean_stride_period(j)*frame_ms total_path_length(j)*scale/512 landing_mean(j,1)*scale/512 landing_mean(j,2)*scale/512 takeoff_mean(j,1)*scale/512 takeoff_mean(j,2)*scale/512 landing_std(j,1)*scale/512 landing_std(j,2)*scale/512 takeoff_std(j,1)*scale/512 takeoff_std(j,2)*scale/512 area(j)*scale/512*scale/512 area_over_path(j)*scale/512 domain_length(j)*scale/512 domain_width(j)*scale/512 footprint_dev(j)]);
    % fprintf(fileID,'\t %0.3f', [100*mov_percentage(j) mean_stride_period(j)*frame_ms total_path_length(j)*scale/512 landing_mean(j,1)*scale/512 landing_mean(j,2)*scale/512 takeoff_mean(j,1)*scale/512 takeoff_mean(j,2)*scale/512 l_std(j)*scale/512 t_std(j)*scale/512 area(j)*scale/512*scale/512 area_over_path(j)*scale/512 domain_length(j)*scale/512 domain_width(j)*scale/512 footprint_dev(j)]);
    fprintf(fileID,'\n');
end
fclose(fileID);
axis([-150 150 -150 150]);
axis square
title('Leg trajectories in body-centered frame of reference');

export_fig([output_dir 'LegDomain.pdf'], '-pdf','-transparent');

fileID = fopen([output_dir 'LegDomainOverlap.xlsx'],'w');
for j = 1:nlegs
     if(missing(j) || nstrides(j)==0)
        continue;
    end
    for k = j+1:nlegs
        if(missing(k) || nstrides(k)==0)
            continue;
        end
        P1.x = nx(convex{j}+startframe-1,j);
        P1.y = ny(convex{j}+startframe-1,j);
        P1.hole = 0;
        P2.x = nx(convex{k}+startframe-1,k);
        P2.y = ny(convex{k}+startframe-1,k);
        P2.hole = 0;
        if(isempty(intersections(P1.x,P1.y,P2.x,P2.y)))
            overlap = 0;
        else
            % Returns the intesersection polygon between the convex polygons for leg j and k
            C = PolygonClip(P1,P2,1);
            overlap = polyarea(C.x,C.y) * (scale/512)^2;
        end
         fprintf(fileID,'%s \t %0.3f \n', ['Overlap ' legs_id{j} '&' legs_id{k}], overlap);
    end
end
fclose(fileID);


%% Estimate leg lengths
% imroi = imread([seg_dir 'roi_' num2str(1) '.png']);
% imsize = size(imroi);
% leglengths = zeros(endframe-startframe+1,nlegs);

% 
% for i = startframe:endframe;
%     imroi = imread([seg_dir 'roi_' num2str(i) '.png']);
%     imseg = imread([seg_dir 'img_' num2str(i) '.png']);
%     imlegs = imseg(:,:,1)==255;
%     imbody = imroi-imlegs;
%     legs = bwconncomp(imlegs);
%     idxlegs = zeros(nlegs,2);
%     idxconncomp = zeros(nlegs,1);
%     for b = 1 :  legs.NumObjects
%         imglegs{b} = zeros(imsize);
%         imglegs{b}(legs.PixelIdxList{b})=1;
%     end
%     for j = 1:nlegs
%         overallMinDistance = inf;
%         for b = 1 : legs.NumObjects     
%            [legY, legX] = ind2sub(imsize, legs.PixelIdxList{b});
%            distances = sqrt((legX - x(i,j)).^2 + (legY - y(i,j)).^2);
%            [minDistance, indexOfMin] = min(distances);
%            if minDistance < overallMinDistance
%                overallMinDistance = minDistance;
%                idxlegs(j,:) = [legX(indexOfMin), legY(indexOfMin)];
%                idxconncomp(j) = b;
%            end
%         end
%         D = bwdistgeodesic(imbody+imglegs{idxconncomp(j)}>0,[idxlegs(j,1)], [idxlegs(j,2)],'quasi-euclidean');
%         leglengths(i-startframe+1,j) = D(round(CoM(i,2)+0.2*body_length),round(CoM(i,1)));
%     end
%     
%     % geodesic distance between middle legs
%     D = bwdistgeodesic(imbody+imglegs{idxconncomp(2)}+imglegs{idxconncomp(5)}>0,[idxlegs(2,1)], [idxlegs(2,2)],'quasi-euclidean');
%     middlelegs_dist(i-startframe+1) = D(idxlegs(5,2),idxlegs(5,1));
%     % horizontal spread distance across middle legs
%     if (norm_trajectory(i-startframe+1,2,1)~=0 && norm_trajectory(i-startframe+1,5,1)~=0)
%         middlelegs_spread(i-startframe+1) = abs(norm_trajectory(i-startframe+1,2,1) -norm_trajectory(i-startframe+1,5,1));
%     else
%         middlelegs_spread(i-startframe+1) = nan;
%     end
% end
% leglengths2 = leglengths;
% leglengths2(find(isnan(leglengths))) = 0;
% leglengths2(find(isinf(leglengths))) = 0;
% fileID = fopen([output_dir 'Leglengths.xlsx'],'w');
% for j = 1: nlegs
%     fprintf(fileID,'\t%s',legs_id{j});
% end
% fprintf(fileID,'\n');
% fprintf(fileID, '%s\t', 'Mean (mm)');
% fclose(fileID);
% mean_ll = zeros (1, nlegs);
% std_ll = zeros (1, nlegs);
% for j = 1 : nlegs
%     A = leglengths2(:,j);
%     mean_ll(j) = mean(A(A>0));
%     std_ll(j) = std(A(A>0));
% end
% dlmwrite([output_dir 'Leglengths.xlsx'],mean_ll*scale/512,'delimiter','\t','-append');
% fileID = fopen([output_dir 'Leglengths.xlsx'],'a');
% fprintf(fileID,'\n');
% fprintf(fileID, '%s\t', 'Std (mm)');
% fclose(fileID);
% dlmwrite([output_dir 'Leglengths.xlsx'],std_ll*scale/512,'delimiter','\t','-append');
% fileID = fopen([output_dir 'Leglengths.xlsx'],'a');
% fprintf(fileID,'\n');
% fprintf(fileID,'Frame #\n');
% fclose(fileID);
% dlmwrite([output_dir 'Leglengths.xlsx'],[(startframe:endframe)', leglengths*scale/512],'delimiter','\t','-append');
% 
% middlelegsAEP = sqrt((landing_mean(2,1)-landing_mean(5,1))^2+(landing_mean(2,2)-landing_mean(5,2))^2);
% middlelegsPEP = sqrt((takeoff_mean(2,1)-takeoff_mean(5,1))^2+(takeoff_mean(2,2)-takeoff_mean(5,2))^2);
% 
% middlelegs_dist2 = middlelegs_dist;
% middlelegs_spread2 = middlelegs_spread;
% middlelegs_dist2(find(isnan(middlelegs_dist))) = 0;
% middlelegs_dist2(find(isinf(middlelegs_dist))) = 0;
% middlelegs_spread2(find(isnan(middlelegs_spread))) = 0;
% middlelegs_spread2(find(isnan(middlelegs_spread))) = 0;
% 
% fileID = fopen([output_dir 'MiddleLegsSpread.xlsx'],'w');
% fprintf(fileID, '%s\t%s\t%s\t%s\t\n', 'Mean Shortest Path Distance (mm)', 'Mean Horizontal Spread (mm)', 'Mean AEP Distance (mm)', 'Mean PEP Distance (mm)');
% fprintf(fileID, '%0.3f\t', mean(middlelegs_dist2(middlelegs_dist2>0)) * scale/512, mean(middlelegs_spread2(middlelegs_spread2>0)) * scale / 512, middlelegsAEP * scale/512, middlelegsPEP * scale/512);
% fprintf(fileID,'\n\n');
% fprintf(fileID, '%s\t', 'Frame #', 'Shortest Path Distance (mm)', 'Horizontal Spread (mm)');
% fprintf(fileID,'\n');
% fclose(fileID);
% dlmwrite([output_dir 'MiddleLegsSpread.xlsx'],[(startframe:endframe)', (middlelegs_dist*scale/512)', (middlelegs_spread*scale/512)'],'delimiter','\t','-append');

middlelegsAEP = sqrt((landing_mean(2,1)-landing_mean(5,1))^2+(landing_mean(2,2)-landing_mean(5,2))^2);
middlelegsPEP = sqrt((takeoff_mean(2,1)-takeoff_mean(5,1))^2+(takeoff_mean(2,2)-takeoff_mean(5,2))^2);
stancewidth = (middlelegsAEP + middlelegsPEP)/2;

fileID = fopen([output_dir 'StanceWidth.xlsx'],'w');
fprintf(fileID, '%s\t%s\t%s\n', 'Mean AEP Distance (mm)', 'Mean PEP Distance (mm)', 'Stance Width (mm)');
fprintf(fileID, '%0.3f\t', middlelegsAEP * scale/512, middlelegsPEP * scale/512, stancewidth * scale / 512);
fclose(fileID);


%% Estimate angle subtended by leg from vertical axis
% legangles = zeros(endframe-startframe+1,nlegs);
% 
% for i = startframe:endframe
%     for j = 1:3
%         if (nx(i,j)==0)
%             legangles(i-startframe+1,j) = nan;
%         else
%             legangles(i-startframe+1,j) = -mod(-atan2d(nx(i,j),ny(i,j)-0.157*bodylength)+360,360);
%         end
%     end
%     for j = 4:6
%         if (nx(i,j)==0)
%             legangles(i-startframe+1,j) = nan;
%         else
%             legangles(i-startframe+1,j) = mod(atan2d(nx(i,j),ny(i,j)-0.157*bodylength)+360,360);
%         end
%     end
% end
% 
% if ~exist([output_dir 'LegAngles'], 'dir')
%    mkdir([output_dir 'LegAngles']);
% end
% 
% fileID = fopen([output_dir 'LegAngles/LegAngles.xlsx'],'w');
% fprintf(fileID,'%s \t %s \t', 'Time(ms)', legs_id{:});
% fprintf(fileID,'\n');
% fclose(fileID);
% dlmwrite([output_dir 'LegAngles/LegAngles.xlsx'],[[1:nframes]'*frame_ms, legangles],'delimiter','\t','-append', 'precision','%.1f');
% 
% PEP_angles = cell(6,1);
% AEP_angles = cell(6,1);
% angle_variations = cell(6,1);
% fileID = fopen([output_dir 'LegAngles/StrideVariation.xlsx'],'w');
% for j = 1:size(stride_start,1)
%     fprintf(fileID,'%s \n', legs_id{j});
%     fprintf(fileID,'%s \t', 'Stride #', 'PEP Angle', 'AEP Angle', 'Angle Swiped');
%     fprintf(fileID,'\n');
%     for i = 1:size(stride_start,2)
%         if (stride_start(j,i)==0)
%             continue;
%         end
%         PEP_angles{j}(i) = legangles(stride_start(j,i),j);
%         AEP_angles{j}(i) = legangles(stride_end(j,i),j);
%         angle_variations{j}(i) = abs(legangles(stride_end(j,i),j)-legangles(stride_start(j,i),j));
%         fprintf(fileID,'%d \t', i, PEP_angles{j}(i), AEP_angles{j}(i), angle_variations{j}(i));
%         fprintf(fileID,'\n');
%     end
% end
% fprintf(fileID,'\n');
% fprintf(fileID,'\n');
% 
% fprintf(fileID,'%s \t', 'Leg', 'Mean PEP Angle', 'Std PEP Angle', 'Mean AEP Angle', 'Std AEP Angle', 'Mean Angle Swiped', 'Std Angle Swiped');
% fprintf(fileID,'\n');
% for j = 1:nlegs
%     fprintf(fileID,'%s \t %d \t', legs_id{j}, mean(PEP_angles{j}), std(PEP_angles{j}), mean(AEP_angles{j}), std(AEP_angles{j}), mean(angle_variations{j}), std(angle_variations{j}));
%     fprintf(fileID,'\n');
% end
% fclose(fileID);


%% Project onto PCA eigenvectors to compute the domain length/width
function [domain_length, domain_width] = projectpca(points)
    m = mean(points);
    v = pca(points);
    n = size(points,1);
    p_v1 = zeros(n,1);
    p_v2 = zeros(n,1);
    for i = 1 : n
        p_v1(i) = projectaxis(points(i,:), v(:,1), m);
        p_v2(i) = projectaxis(points(i,:), v(:,2), m);
    end
    domain_length = max(p_v1) - min(p_v1);
    domain_width = max(p_v2) - min(p_v2);
    
    
function p_proj = projectaxis(p, v, m, axis)
    p_proj = p*v-m*v;