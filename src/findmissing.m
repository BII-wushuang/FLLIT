function findmissing()
data_dir = uipickfiles('FilterSpec',[pwd '/Results/Tracking']);
fileID = fopen('missing_tips.csv','w');
fprintf(fileID,'%s \t %s \t %s \t %s \t %s \t %s \t %s \t %s', 'Dataset', '#frames', 'L1', 'L2', 'L3', 'L4', 'L5', 'L6');
fprintf(fileID,'\n');
fclose(fileID);

for i = 1 : length(data_dir)
    ProcessFolder(data_dir{i});
end
end

function ProcessFolder(data_dir)
    pos_bs = strfind(data_dir,'Tracking');
    sub_dir = data_dir(pos_bs(end)+length('Tracking'):length(data_dir));
    if(exist(data_dir) && isempty(dir([data_dir '/trajectory.mat'])) )
        dir_fn = dir(data_dir);
        for i_dir = 1:length(dir_fn)-2
            ProcessFolder([data_dir '/' dir_fn(i_dir+2).name]);
            pause(1);
        end
    else
        data_dir = ['./Results/Tracking' sub_dir];
        data_name_pos = strfind(data_dir,'/');
        data_name = data_dir(data_name_pos(end)+1:end);
        load([data_dir '/trajectory.mat']);
        output{1} = data_name;
        output{2} = length(trajectory);
        for j =1:6
            output{j+2} = length(find(trajectory(:,j,1)==0));
        end
        fileID = fopen('missing_tips.csv','a');
        formatSpec = '%s \t %d \t %d \t %d \t %d \t %d \t %d \t %d \t %d';
        fprintf(fileID,formatSpec,output{:});
        fprintf(fileID,'\n');
        fclose(fileID);
    end
end