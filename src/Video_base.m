%Plot the trajectory of the fly's legs
function Video(data_dir,fps,skip,startframe,endframe)

pos_bs = strfind(data_dir,'Data');
sub_dir = data_dir(pos_bs(end)+length('Data'):length(data_dir));
data_dir = ['./Data' sub_dir '/'];
if(~isempty(dir([data_dir '*.tif'])))
    img_list = dir([data_dir '*.tif']);
else
    img_list = dir([data_dir '*.bmp']);
end

load(['./Results/Tracking/' sub_dir '/trajectory.mat']);
load(['./Results/Tracking/' sub_dir '/norm_trajectory.mat']);
load(['./Results/Tracking/' sub_dir '/CoM.mat']);

myVideo = VideoWriter(['./Results/Tracking' sub_dir '/Video.avi']);
myVideo.FrameRate = fps;
myVideo.Quality = 100;
nlegs = size(trajectory,2);
if (nlegs == 6)
    legs_id = {'L1', 'L2', 'L3', 'R1', 'R2', 'R3'};
    colors = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 0 1 1; 1 0 1];
else
    legs_id = {'L1', 'L2', 'L3', 'L4', 'R1', 'R2', 'R3', 'R4'};
    trajectory = trajectory(:,[1 2 3 7 4 5 6 8],:);
    norm_trajectory = norm_trajectory(:,[1 2 3 7 4 5 6 8],:);
    colors = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 0 1 1; 1 0 1; 1 0.5 0; 1 0 0.5];
end

start = zeros(nlegs,1);
for j = 1:nlegs;
    start(j) = find(trajectory(:,j,1),1);
end
startf = max(start);
endf = length(CoM);

for j =1:size(trajectory,2)
    nonzero = find(norm_trajectory(:,j,1));
    for i = startf : endf
        if (trajectory(i,j,1) ~= 0)
            x(i,j) = trajectory(i,j,1);
            y(i,j) = trajectory(i,j,2);
            nx(i,j) = norm_trajectory(i,j,1);
            ny(i,j) = norm_trajectory(i,j,2);
        else
            previdx = nonzero(find(nonzero<i, 1));
            nextidx = nonzero(find(nonzero>i, 1));
            if(isempty(previdx))
                x(i,j) = trajectory(nextidx,j,1);
                y(i,j) = trajectory(nextidx,j,2);
                nx(i,j) = norm_trajectory(nextidx,j,1);
                ny(i,j) = norm_trajectory(nextidx,j,2);
            elseif (~isempty(nextidx))
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

for i = startf : endf
    if (i == startf)
        dist(i) = 0;
    else
        dist(i) = dist(i-1) + (CoM(i,1) - CoM(i-1,1))*sind(CoM(i-1,3)) + (CoM(i,2) - CoM(i-1,2))*-cosd(CoM(i-1,3));
    end
end

if(isempty(strfind(sub_dir, '/')))
    pos_bs = strfind(sub_dir, '\');
else
    pos_bs = strfind(sub_dir, '/');
end
Vtitle = sub_dir(pos_bs(end)+1:end);

fig = figure(1);
fig.Position = [0 0 1920 1080];
fig.Resize = 'on';
traj(1) = getframe(fig);

counter = 1;
for i = startframe : skip : endframe
    counter = counter+1;
    if mod(i,50) == 0
        fprintf('Processing frame %d\n',i);
    end
    
    t = startframe : skip : i;
    nt = t;
    
    figure(1);
    img = imread([data_dir img_list(i).name]);
    hold on
    subplot(4,4,[2:3 6:7])
    imshow(img);
    axis square
    hold off
    
    img_norm = imtranslate(img, [255 - CoM(i,1), 255 - CoM(i,2)]);
    img_norm = imrotate(img_norm, CoM(i,3));
    img_norm = imcrop(img_norm, [size(img_norm,1)/2-150 size(img_norm,2)/2-150 300 300]);
    hold on
    subplot(4,4,[10:11 14:15])
    imshow(img_norm);
    axis square
    hold off
    
    hold on
    % Superposed trajectory with fly
    subplot(4,4,[2:3 6:7])
    for j = 1: size(trajectory,2)
        p(j) = plot (x(t,j), y(t,j));
    end
    scatter (CoM(i,1),CoM(i,2),'w*');
    
    title(Vtitle,'Interpreter', 'none');
    axis([0 512 0 512]);
    axis square
    set(gca,'Ydir','reverse')
    
    hold on
    % Superposed trajectory with fly in normalised plane
    subplot(4,4,[10:11 14:15])
    for j = 1: size(trajectory,2)
        np(j) = plot (150+nx(t,j), 150+ny(t,j));
    end
    title('Leg motion in normalised plane');
    axis square
    
    % CoM trajectory in x-y plane
    subplot(4,4,[1 5])
    plot (CoM(t,1), 512-CoM(t,2));
    title('CoM trajectory in x-y plane');
    axis([0 512 0 512]);
    axis square
    
    % Forward displacement of fly
    subplot(4,4,[4 8])
    plot(nt, dist(t));
    title('Forward displacement (direction it is heading)');
    xlabel('time (ms)');
    ylabel('displacement (pixels)');
    axis square
    
    % Vertical displacement of legs
    subplot(4,4,[9 13])
    npy = plot (nt, -ny(t,1:size(trajectory,2)));
    title('Vertical displacement of legs');
    xlabel('time (ms)');
    ylabel('displacement (pixels)');
    axis square
    set(gca,'Color',[0.5 0.5 0.5]);
    
    % Lateral displacement of legs
    subplot(4,4,[12 16])
    npx = plot (nt, nx(t,1: size(trajectory,2)));
    title('Lateral displacement of legs');
    xlabel('time (ms)');
    ylabel('displacement (pixels)');
    axis square
    set(gca,'Color',[0.5 0.5 0.5]);
    
    hold off
     
    for j = 1 : size(trajectory,2)
        p(j).Color = [colors(j,1) colors(j,2) colors(j,3)];
        p(j).LineStyle = '--';
        np(j).Color = [colors(j,1) colors(j,2) colors(j,3)];
        np(j).LineStyle = '--';
        npx(j).Color = [colors(j,1) colors(j,2) colors(j,3)];
        npy(j).Color = [colors(j,1) colors(j,2) colors(j,3)];
    end
    
    drawnow
    traj(counter) = getframe(fig);
end

open(myVideo);
writeVideo(myVideo, traj);
close(myVideo);