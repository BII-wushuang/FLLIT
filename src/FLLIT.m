% FLLIT program: GUI incorporating the program workflow
% This GUI is created with the MATLAB GUIDE feature.
% The pushbuttons call functions to perform individual tasks such as
% segmentation, tracking or data processing.
function varargout = FLLIT(varargin)
% FLLIT MATLAB code for FLLIT.fig
%      FLLIT, by itself, creates a new FLLIT or raises the existing
%      singleton*.
%
%      H = FLLIT returns the handle to a new FLLIT or the handle to
%      the existing singleton*.
%
%      FLLIT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FLLIT.M with the given input arguments.
%
%      FLLIT('Property','Value',...) creates a new FLLIT or raises
%      the existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FLLIT_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FLLIT_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FLLIT

% Last Modified by GUIDE v2.5 11-Sep-2018 11:33:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FLLIT_OpeningFcn, ...
                   'gui_OutputFcn',  @FLLIT_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});

end
% End initialization code - DO NOT EDIT


%% --- Executes just before FLLIT is made visible.
function FLLIT_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   Command line arguments to FLLIT (see VARARGIN)

% Addpath
try
    addpath(genpath('./KernelBoost-v0.1/'));
catch
end

% Choose default Command line output for FLLIT
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
% 'gcf' refers to 'get current figure' and in this context all the updated
% data are feed into the 'hMainGui'. The 'hMainGui' can then provide the
% data as global variables to the other functions that are called in this GUI
setappdata(0, 'hMainGui', gcf);

S.fh = handles.figMain;
% Fix the central figure axis in pixel coordinates
S.ax = gca;
set(S.ax,'unit','pix','position',[(handles.figMain.Position(3)-512)/2 handles.figMain.Position(4)-512 512 512]);
S.XLM = get(S.ax,'xlim');
S.YLM = get(S.ax,'ylim');
S.AXP = get(S.ax,'pos');
S.DFX = diff(S.XLM);
S.DFY = diff(S.YLM);
S.tx = handles.Position_Text;
% The fh_wbmcfn function captures motion of the mouse cursor and the cursor
% location is displayed in the GUI
set(S.fh,'windowbuttonmotionfcn',{@fh_wbmfcn,S}) 

% The clicker function captures double clicks of the mouse and is used for
% annotating / labelling of the leg tip locations
set(handles.figMain, 'WindowButtonDownFcn', {@clicker,handles});

initialize_gui(hObject, handles, false);


%% --- Outputs from this function are returned to the Command line.
function varargout = FLLIT_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default Command line output from handles structure
varargout{1} = handles.output;


%% --------------------------------------------------------------------
function initialize_gui(fig_handle, handles, isreset)

% Update handles structure
guidata(handles.figMain, handles);


%% --- Executes on button press in Adjust_Prediction.
function Adjust_Prediction_Callback(hObject, eventdata, handles)
% hObject    handle to Adjust_Prediction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% First, obtain the global variables stored in 'hMainGui': named the
% current frame in the tracking progess as well as currently stored
% trajectory data of the legs
hMainGui = getappdata(0, 'hMainGui');
i = getappdata(hMainGui, 'current_frame');
trajectory = getappdata(hMainGui, 'trajectory');

nLegs = 8 - handles.L4.Value - handles.R4.Value;
if (nLegs > 6)
    colors = [1 0 0; 0 1 0; 0 0 1; 1 0 0.5; 1 1 0; 0 1 1; 1 0 1; 1 0.5 0];
    legs_id = {'L1', 'L2', 'L3','L4', 'R1', 'R2', 'R3','R4'};
else
    colors = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 0 1 1; 1 0 1; 1 0 0.5; 1 0.5 0];
    legs_id = {'L1', 'L2', 'L3', 'R1', 'R2', 'R3'};
end

if (isempty(i))
    % If current frame is empty, it means that tracking has not been run 
    handles.Text_Display.String = 'Adjustments are to correct tracking errors, please run tracking first.';
elseif (strcmp(hObject.String, 'Adjust Prediction'))
    % Clicking on the adjust prediction button will result in itself being
    % displayed as the 'Exit without Saving' button. Once the adjust
    % prediction mode is activated, double clicking on the image will be
    % registered and return the tip location of the selected leg. 
    hObject.String = 'Exit without saving';
    % A new button 'Save and Exit' will also appear. The handle to this
    % 'Save and Exit' button is handles.Annotation.
    handles.Annotation.Enable = 'On';
    handles.Annotation.Visible = 'On';
    handles.Annotation.String = 'Save and Exit';
    handles.ConsoleUpdate.Enable = 'On';
    handles.ConsoleUpdate.Visible = 'On';
    handles.prevframe.Enable = 'Off';
    handles.prevframe.Visible = 'Off';
    handles.nextframe.Enable = 'Off';
    handles.nextframe.Visible = 'Off';
    if (~strcmp(handles.Tracking.String, 'Initial'))
        handles.Tracking.String = 'Resume';
    end

    set(handles.Text_Display,'String',['Entering adjustment mode: select a leg and double click on figure to assign tip location.']);
elseif (strcmp(hObject.String, 'Exit without saving'))
    % if user decides to exit without saving, all user input on the current frame will be
    % discarded
    hObject.String = 'Adjust Prediction';
    handles.Annotation.String = 'Annotate';
    handles.Annotation.Enable = 'Off';
    handles.Annotation.Visible = 'Off';
    handles.ConsoleUpdate.Enable = 'Off';
    handles.ConsoleUpdate.Visible = 'Off';
    handles.prevframe.Enable = 'On';
    handles.prevframe.Visible = 'On';
    handles.nextframe.Enable = 'On';
    handles.nextframe.Visible = 'On';
    handles.Tracking_Text.String = '';
    for kk = 1:nLegs
        if(trajectory(i,kk,1)~=0)
            handles.Tracking_Text.String = [handles.Tracking_Text.String; sprintf('Leg %2s: %3d %3d', legs_id{kk}, trajectory(i,kk,1), trajectory(i,kk,2))]; 
        else
            handles.Tracking_Text.String = [handles.Tracking_Text.String; sprintf('Leg %2s:        ', legs_id{kk})];
            for j = 1: 2 : length(handles.figMain.CurrentAxes.Children)-2
               if(strcmp(handles.figMain.CurrentAxes.Children(j).Type,'text') && strcmp(handles.figMain.CurrentAxes.Children(j).String,legs_id{kk}))
                   delete(handles.figMain.CurrentAxes.Children(j:j+1));
               end
            end
        end
    end
    for i = 1 : nLegs
        x = str2num(handles.Tracking_Text.String(i,9:11));
        y = str2num(handles.Tracking_Text.String(i,13:15));
        for j = 1: 2 : length(handles.figMain.CurrentAxes.Children)-2
           if(strcmp(handles.figMain.CurrentAxes.Children(j).Type,'text') && strcmp(handles.figMain.CurrentAxes.Children(j).String,legs_id{i}))
               delete(handles.figMain.CurrentAxes.Children(j:j+1));
               if(x~=0)
               hold on
               scatter(x,y,'w');
               text(x,y,legs_id{i},'Color',[colors(i,1) colors(i,2) colors(i,3)],'FontSize',14);
               hold off
               end
           end
        end
    end
end
handles.ManualInitiate.Value = 0;


%% --- Executes on button press in ConsoleUpdate.
function ConsoleUpdate_Callback(hObject, eventdata, handles)
% hObject    handle to ConsoleUpdate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hMainGui = getappdata(0, 'hMainGui');

nLegs = 8 - handles.L4.Value - handles.R4.Value;
if (nLegs > 6)
    colors = [1 0 0; 0 1 0; 0 0 1; 1 0 0.5; 1 1 0; 0 1 1; 1 0 1; 1 0.5 0];
    legs_id = {'L1', 'L2', 'L3','L4', 'R1', 'R2', 'R3','R4'};
else
    colors = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 0 1 1; 1 0 1; 1 0 0.5; 1 0.5 0];
    legs_id = {'L1', 'L2', 'L3', 'R1', 'R2', 'R3'};
end

for i = 1 : nLegs
    x = str2num(handles.Tracking_Text.String(i,9:11));
    y = str2num(handles.Tracking_Text.String(i,13:15));
    for j = 1 : length(handles.figMain.CurrentAxes.Children)-2
       if(strcmp(handles.figMain.CurrentAxes.Children(j).Type,'text') && strcmp(handles.figMain.CurrentAxes.Children(j).String,legs_id{i}))
           delete(handles.figMain.CurrentAxes.Children(j:j+1));
           hold on
           scatter(x,y,'w');
           text(x,y,legs_id{i},'Color',[colors(i,1) colors(i,2) colors(i,3)],'FontSize',14);
           hold off
       end
    end
end


%% --- Executes on button press in Annotation.
function Annotation_Callback(hObject, eventdata, handles)
% hObject    handle to Annotation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% This allows the user input to relabel the legs on the current frame. This
% will then update the data in the trajectory as well as the normalised
% trajectory.

hMainGui = getappdata(0, 'hMainGui');
data_dir = getappdata(hMainGui, 'data_dir');
pos_bs = strfind(data_dir,'Data');
sub_dir = data_dir(pos_bs(end)+length('Data'):length(data_dir));
output_dir = ['./Results/Tracking/' sub_dir '/'];
if(~exist(output_dir))
    mkdir(output_dir);
end
CoM = getappdata(hMainGui, 'CoM');
trajectory = getappdata(hMainGui, 'trajectory');
norm_trajectory = getappdata(hMainGui, 'norm_trajectory');
current_frame = getappdata(hMainGui,'current_frame');

nLegs = 8 - handles.L4.Value - handles.R4.Value;

trajectory(current_frame,1,1) = str2num(handles.Tracking_Text.String(1,9:11));
trajectory(current_frame,1,2) = str2num(handles.Tracking_Text.String(1,13:15));
norm_trajectory(current_frame,1,1) = cosd(CoM(current_frame,3))*((trajectory(current_frame,1,1)-CoM(current_frame,1))) + sind(CoM(current_frame,3))*((trajectory(current_frame,1,2)-CoM(current_frame,2)));

if(strcmp(handles.Tracking.String, 'Initial') && norm_trajectory(current_frame,1,1)>0)
        CoM(current_frame,3) = mod(CoM(current_frame,3)+180,360);
        setappdata(hMainGui,'CoM',CoM);
        % find body_length
        i = current_frame;
        img = imread([seg_dir 'roi_' num2str(i) '.png']);
        img_norm = imtranslate(img, [255 - CoM(i,1), 255 - CoM(i,2)]);
        img_norm = imrotate(img_norm, CoM(i,3));
        img_norm = imcrop(img_norm, [size(img_norm,1)/2-150 size(img_norm,2)/2-150 300 300]);
        [Y,X] = find(img_norm);
        bodylength = max(Y(find(X==150)))-min(Y(find(X==150)));
        handles.BLText.Visible = 'On';
        handles.BLText.Visible = 'On';
        handles.BodyLength.String = num2str(round(bodylength));
        ResetImageSize_Callback(hObject, eventdata, handles);
end

for i = 1: nLegs
    if (~isempty(str2num(handles.Tracking_Text.String(i,9:11))))
       trajectory(current_frame,i,1) = str2num(handles.Tracking_Text.String(i,9:11));
       trajectory(current_frame,i,2) = str2num(handles.Tracking_Text.String(i,13:15));
       norm_trajectory(current_frame,i,1) = cosd(CoM(current_frame,3))*((trajectory(current_frame,i,1)-CoM(current_frame,1))) + sind(CoM(current_frame,3))*((trajectory(current_frame,i,2)-CoM(current_frame,2)));
       norm_trajectory(current_frame,i,2) = -sind(CoM(current_frame,3))*((trajectory(current_frame,i,1)-CoM(current_frame,1))) + cosd(CoM(current_frame,3))*((trajectory(current_frame,i,2)-CoM(current_frame,2)));
    else
       trajectory(current_frame,i,:) = 0;
       norm_trajectory(current_frame,i,:) = 0;
    end
end

setappdata(hMainGui,'trajectory',trajectory);
setappdata(hMainGui,'norm_trajectory',norm_trajectory);
save([output_dir 'CoM.mat'],'CoM');
save([output_dir 'trajectory.mat'],'trajectory');
save([output_dir 'norm_trajectory.mat'],'norm_trajectory');
csvwrite([output_dir 'CoM.csv'],CoM);
csvwrite([output_dir 'trajectory.csv'],trajectory);
csvwrite([output_dir 'norm_trajectory.csv'],norm_trajectory);

hObject.String = 'Annotate';
hObject.Enable = 'Off';
hObject.Visible = 'Off';
handles.ConsoleUpdate.Enable = 'Off';
handles.ConsoleUpdate.Visible = 'Off';
handles.prevframe.Enable = 'On';
handles.prevframe.Visible = 'On';
handles.nextframe.Enable = 'On';
handles.nextframe.Visible = 'On';
handles.Adjust_Prediction.String = 'Adjust Prediction';
handles.Tracking.String = 'Resume';


%% --- Executes on button press in prevframe.
function prevframe_Callback(hObject, eventdata, handles)
% hObject    handle to prevframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hMainGui = getappdata(0, 'hMainGui');
data_dir = getappdata(hMainGui, 'data_dir');
pos_bs = strfind(data_dir,'Data');
sub_dir = data_dir(pos_bs(end)+length('Data'):length(data_dir));
seg_dir = ['./Results/SegmentedImages' sub_dir '/'];
data_dir = ['./Data' sub_dir '/'];
if(~isempty(dir([data_dir '*.tif'])))
    img_list = dir([data_dir '*.tif']);
else
    img_list = dir([data_dir '*.bmp']);
end

setappdata(hMainGui,'current_frame', getappdata(hMainGui, 'current_frame') - 1);
i = getappdata(hMainGui, 'current_frame');
handles.text.String = ['Current frame: ' num2str(i)];
handles.GoToDropdown.Value = i;

display = handles.ChooseDisplay.Value;
switch display
    case 1
        handles.GoToDropdown.String = 1:length(img_list);
        imshow([data_dir img_list(i).name]);
    case 2
        seg_list = dir([seg_dir 'img_*.png']);
        handles.GoToDropdown.String = 1:length(seg_list);
        imshow([seg_dir 'img_' num2str(i) '.png']);
    case 3
        handles.Tracking.String = 'Resume';
        trajectory = getappdata(hMainGui, 'trajectory');

        nLegs = 8 - handles.L4.Value - handles.R4.Value;
        if (nLegs > 6)
            colors = [1 0 0; 0 1 0; 0 0 1; 1 0 0.5; 1 1 0; 0 1 1; 1 0 1; 1 0.5 0];
            legs_id = {'L1', 'L2', 'L3','L4', 'R1', 'R2', 'R3','R4'};
        else
            colors = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 0 1 1; 1 0 1; 1 0 0.5; 1 0.5 0];
            legs_id = {'L1', 'L2', 'L3', 'R1', 'R2', 'R3'};
        end

        if(~isempty(trajectory))
            handles.text.String = ['Current frame: ' num2str(i)];

            handles.Tracking_Text.String = '';
            imshow([data_dir img_list(i).name]);
            hold on;
            for kk = 1:nLegs
                if(eval(['handles.' legs_id{kk} '.Value']))
                    handles.Tracking_Text.String = [handles.Tracking_Text.String; sprintf('Leg %2s:        ', legs_id{kk})];
                    continue;
                end
                if(trajectory(i,kk,1)~=0)
                    handles.Tracking_Text.String = [handles.Tracking_Text.String; sprintf('Leg %2s: %3d %3d', legs_id{kk}, trajectory(i,kk,1), trajectory(i,kk,2))]; 
                    scatter(trajectory(i,kk,1),trajectory(i,kk,2),'w');
                    text(trajectory(i,kk,1),trajectory(i,kk,2),legs_id{kk},'Color',[colors(kk,1) colors(kk,2) colors(kk,3)],'FontSize',14);
                else
                    handles.Tracking_Text.String = [handles.Tracking_Text.String; sprintf('Leg %2s:        ', legs_id{kk})];
                end
            end
            hold off;
            drawnow
            
            % Update Display
            nLegs = 8 - handles.L1.Value - handles.L2.Value - handles.L3.Value - handles.R1.Value - handles.R2.Value - handles.R3.Value - handles.L4.Value - handles.R4.Value;
            if(length(str2num(handles.Tracking_Text.String(:,9:11))) < nLegs)
                set(handles.Missing, 'BackgroundColor', [1 0 0]); 
            else
                set(handles.Missing, 'BackgroundColor', [0 1 0]); 
            end
        end
end


%% --- Executes on button press in nextframe.
function nextframe_Callback(hObject, eventdata, handles)
% hObject    handle to nextframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hMainGui = getappdata(0, 'hMainGui');
data_dir = getappdata(hMainGui, 'data_dir');
pos_bs = strfind(data_dir,'Data');
sub_dir = data_dir(pos_bs(end)+length('Data'):length(data_dir));
seg_dir = ['./Results/SegmentedImages' sub_dir '/'];
data_dir = ['./Data' sub_dir '/'];
if(~isempty(dir([data_dir '*.tif'])))
    img_list = dir([data_dir '*.tif']);
else
    img_list = dir([data_dir '*.bmp']);
end

setappdata(hMainGui,'current_frame', getappdata(hMainGui, 'current_frame') + 1);
i = getappdata(hMainGui, 'current_frame');
handles.text.String = ['Current frame: ' num2str(i)];
handles.GoToDropdown.Value = i;

display = handles.ChooseDisplay.Value;
switch display
    case 1
        handles.GoToDropdown.String = 1:length(img_list);
        imshow([data_dir img_list(i).name]);
    case 2
        seg_list = dir([seg_dir 'img_*.png']);
        handles.GoToDropdown.String = 1:length(seg_list);
        imshow([seg_dir 'img_' num2str(i) '.png']);
    case 3
        CoM = getappdata(hMainGui, 'CoM');
        trajectory = getappdata(hMainGui, 'trajectory');
        norm_trajectory = getappdata(hMainGui, 'norm_trajectory');
        se = strel('disk',4);
        skip = 1;
        max_distance = 20;
        
        nLegs = 8 - handles.L4.Value - handles.R4.Value;
        if (nLegs > 6)
            colors = [1 0 0; 0 1 0; 0 0 1; 1 0 0.5; 1 1 0; 0 1 1; 1 0 1; 1 0.5 0];
            legs_id = {'L1', 'L2', 'L3','L4', 'R1', 'R2', 'R3','R4'};
        else
            colors = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 0 1 1; 1 0 1; 1 0 0.5; 1 0.5 0];
            legs_id = {'L1', 'L2', 'L3', 'R1', 'R2', 'R3'};
        end
        
        if (handles.ViewerMode.Value)
            if(~isempty(trajectory))
                handles.text.String = ['Current frame: ' num2str(i)];

                handles.Tracking_Text.String = '';
                imshow([data_dir img_list(i).name]);
                hold on;
                for kk = 1:nLegs
                    if(eval(['handles.' legs_id{kk} '.Value']))
                        handles.Tracking_Text.String = [handles.Tracking_Text.String; sprintf('Leg %2s:        ', legs_id{kk})];
                        continue;
                    end
                    if(trajectory(i,kk,1)~=0)
                        handles.Tracking_Text.String = [handles.Tracking_Text.String; sprintf('Leg %2s: %3d %3d', legs_id{kk}, trajectory(i,kk,1), trajectory(i,kk,2))]; 
                        scatter(trajectory(i,kk,1),trajectory(i,kk,2),'w');
                        text(trajectory(i,kk,1),trajectory(i,kk,2),legs_id{kk},'Color',[colors(kk,1) colors(kk,2) colors(kk,3)],'FontSize',14);
                    else
                        handles.Tracking_Text.String = [handles.Tracking_Text.String; sprintf('Leg %2s:        ', legs_id{kk})];
                    end
                end
                hold off;
                drawnow
                
                % Update Display
                nLegs = 8 - handles.L1.Value - handles.L2.Value - handles.L3.Value - handles.R1.Value - handles.R2.Value - handles.R3.Value - handles.L4.Value - handles.R4.Value;
                if(length(str2num(handles.Tracking_Text.String(:,9:11))) < nLegs)
                    set(handles.Missing, 'BackgroundColor', [1 0 0]); 
                else
                    set(handles.Missing, 'BackgroundColor', [0 1 0]); 
                end
            end
        else
            handles.Tracking.String = 'Resume';
            
            Im = (imread([seg_dir 'img_' num2str(i) '.png']));
            Imroi = imread([seg_dir 'roi_' num2str(i) '.png']);

            % Fix image/segmentation irregularities
            Imwork = double(Im(:,:,1) == 255) .* Imroi;
            
            if (nLegs>6)
                Imwork = Imwork - (bwdist(imerode(Imroi,strel('disk',6)))<15);
				Imwork = bwareaopen(Imwork>0,5);
            else
                Imwork = Imwork - (bwdist(imerode(Imroi,se))<10);
				Imwork = bwareaopen(Imwork>0,3);
            end
            
            [locY,locX] = find(imerode(Imroi,se) > 0);

            [centre_of_mass,theta] = findCoM(locX,locY);
            
            if (abs(theta - CoM(i-skip,3)) > 90 && abs(theta - CoM(i-skip,3)) < 270)
                theta = mod(theta + 180,360);
            end

            CoM(i, :) = [centre_of_mass, theta];

            clear raw_tips;
            clear raw_normtips;
            
            [raw_tips, raw_normtips] = findtips(Imwork, Imroi, locX, locY, centre_of_mass, theta);
            
            x_t1 = raw_normtips;
            x_t0 = reshape(norm_trajectory(i-1,:),[nLegs 2]);
            target_indices = hungarianlinker(x_t0, x_t1, max_distance);

            for kk = 1:length(target_indices)
                if (target_indices(kk) > 0)
                    norm_trajectory(i, kk, :) = raw_normtips(target_indices(kk),:);
                    trajectory(i, kk, :) = raw_tips(target_indices(kk),:);
                else
                    norm_trajectory(i, kk, :) = 0;
                    trajectory(i, kk, :) = 0;
                end
            end

            leftoveridx = find(hungarianlinker(raw_normtips, x_t0, max_distance)==-1);
            x_t1 = raw_normtips(leftoveridx,:);
            x_t2 = raw_tips(leftoveridx,:);

            last_seen_tips = reshape(norm_trajectory(i,:),[nLegs 2]);
            unassigned_idx = find(last_seen_tips(:,1) == 0);
            for kk = 1:length(unassigned_idx)
                if(eval(['handles.' legs_id{unassigned_idx(kk)} '.Value']))
                    continue;
                end
                idx = find(norm_trajectory(1:i,unassigned_idx(kk),1),1,'last');
                last_seen_tips(unassigned_idx(kk),:) = reshape(norm_trajectory(idx,unassigned_idx(kk),:), [1 2]);
            end

            if (~isempty(unassigned_idx) && ~isempty(x_t1))
                    C = hungarianlinker(last_seen_tips(unassigned_idx,:), x_t1, 1.25*max_distance);
                    for l = 1:length(C)
                        if (C(l) ~= -1)
                            norm_trajectory(i,unassigned_idx(l), :) = x_t1(C(l), :);
                            trajectory(i,unassigned_idx(l), :) = x_t2(C(l), :);
                        end
                    end
            end

            handles.Tracking_Text.String = '';
            oriIm = im2double(imread([data_dir img_list(i).name]));
            imshow(oriIm);
            hold on;
            for kk = 1:nLegs
                if (eval(['handles.' legs_id{kk} '.Value']))
                        trajectory(i,kk,1) = 0;
                        trajectory(i,kk,2) = 0;
                        norm_trajectory(i,kk,1) = 0;
                        norm_trajectory(i,kk,2) = 0;
                        handles.Tracking_Text.String = [handles.Tracking_Text.String; sprintf('Leg %2s:        ', legs_id{kk})];
                        continue;
                end
                if(trajectory(i,kk,1)~=0)
                    handles.Tracking_Text.String = [handles.Tracking_Text.String; sprintf('Leg %2s: %3d %3d', legs_id{kk}, trajectory(i,kk,1), trajectory(i,kk,2))]; 
                    scatter(trajectory(i,kk,1),trajectory(i,kk,2),'w');
                    text(trajectory(i,kk,1),trajectory(i,kk,2),legs_id{kk},'Color',[colors(kk,1) colors(kk,2) colors(kk,3)],'FontSize',14);
                else
                    handles.Tracking_Text.String = [handles.Tracking_Text.String; sprintf('Leg %2s:        ', legs_id{kk})];
                end
            end
            hold off;
            drawnow
            
            % Update Display
            nLegs = 8 - handles.L1.Value - handles.L2.Value - handles.L3.Value - handles.R1.Value - handles.R2.Value - handles.R3.Value - handles.L4.Value - handles.R4.Value;
            if(length(str2num(handles.Tracking_Text.String(:,9:11))) < nLegs)
                set(handles.Missing, 'BackgroundColor', [1 0 0]); 
            else
                set(handles.Missing, 'BackgroundColor', [0 1 0]); 
            end

            setappdata(hMainGui, 'CoM' , CoM);
            setappdata(hMainGui, 'trajectory' , trajectory);
            setappdata(hMainGui, 'norm_trajectory' , norm_trajectory);
            handles.GoToDropdown.String = 1:find(trajectory,1,'last');
        end
end


%% --- Executes on button press in Select_Folder.
function Select_Folder_Callback(hObject, eventdata, handles)
% hObject    handle to Select_Folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hMainGui = getappdata(0, 'hMainGui');
if(~isempty(getappdata(hMainGui, 'trajectory')))        
    choice = questdlg('Save tracking data?', ...
        'FLLIT', ...
        'Yes', 'No', 'Cancel', 'Cancel');

    % Handle response
    switch choice
        case 'Yes'
                data_dir = getappdata(hMainGui, 'data_dir');
                pos_bs = strfind(data_dir,'Data');
                sub_dir = data_dir(pos_bs(end)+length('Data'):length(data_dir));
                output_dir = ['./Results/SegmentedImages' sub_dir '/'];
                trajectory = getappdata(hMainGui, 'trajectory');
                norm_trajectory = getappdata(hMainGui, 'norm_trajectory');
                CoM = getappdata(hMainGui, 'CoM');
                save([output_dir 'CoM.mat'],'CoM');
                save([output_dir 'trajectory.mat'],'trajectory');
                save([output_dir 'norm_trajectory.mat'],'norm_trajectory');
                csvwrite([output_dir 'CoM.csv'],CoM);
                csvwrite([output_dir 'trajectory.csv'],trajectory);
                csvwrite([output_dir 'norm_trajectory.csv'],norm_trajectory);
            	setappdata(hMainGui, 'trajectory', []);
        case 'No'
            setappdata(hMainGui, 'trajectory', []);
        case 'Cancel'
            return;
    end
end

handles.BLText.Visible = 'Off';
handles.BodyLength.Visible = 'Off';
handles.BodyLength.String = '0';
handles.Segmentation.Enable = 'Off';
handles.Tracking.Enable = 'Off';
handles.Adjust_Prediction.Enable = 'Off';
handles.Video.Enable = 'Off';
handles.DataProcessing.Enable = 'Off';
handles.prevframe.Enable = 'Off';
handles.nextframe.Enable = 'Off';
handles.GoToFrame.Visible = 'Off';
handles.GoToDropdown.Visible = 'Off';
handles.ChooseDisplay.Visible = 'Off';
handles.Tracking_Text.String = '';

data_dir = uigetdir(['./Data']);
pos_bs = strfind(data_dir,'Data');
if(~isempty(pos_bs))
    sub_dir = data_dir(pos_bs(end)+length('Data'):length(data_dir));
    setappdata(hMainGui, 'data_dir', data_dir);
    set(handles.Text_Display, 'String', data_dir);
else
    return;
end
data_dir = ['./Data' sub_dir '/'];
handles.Progress.String = '';
handles.Tracking.String = 'Tracking';

if(~isempty(dir([data_dir '*.tif'])))
    img_list = dir([data_dir '*.tif']);
else
    img_list = dir([data_dir '*.bmp']);
end

if(~isempty(dir([data_dir '*.tif'])) || ~isempty(dir([data_dir '*.bmp'])))
    handles.Segmentation.Enable = 'On';
    handles.Tracking.Enable = 'On';
    handles.Adjust_Prediction.Enable = 'On';
    handles.Video.Enable = 'On';
    handles.DataProcessing.Enable = 'On';
    handles.prevframe.Enable = 'On';
    handles.nextframe.Enable = 'On';
   
    handles.GoToFrame.Visible = 'On';
    handles.GoToDropdown.Visible = 'On';
    handles.ChooseDisplay.Visible = 'On';
    handles.GoToDropdown.String = 1:length(img_list);
    handles.GoToDropdown.Value = 1;
    handles.ChooseDisplay.String = {'Original'};
    handles.ChooseDisplay.Value = 1;
    handles.text.String = 'Current frame: 1';
    setappdata(hMainGui,'current_frame',1);
    handles.Tracking_Text.String = '';
    oriIm = im2double(imread([data_dir img_list(1).name]));
    imshow(oriIm);
else
    handles.Segmentation.Enable = 'Off';
    handles.Tracking.Enable = 'Off';
    handles.Video.Enable = 'Off';
    handles.Adjust_Prediction.Enable = 'Off';
    handles.DataProcessing.Enable = 'Off';
    handles.prevframe.Enable = 'Off';
    handles.nextframe.Enable = 'Off';
    handles.GoToFrame.Visible = 'Off';
    handles.GoToDropdown.Visible = 'Off';
    handles.ChooseDisplay.Visible = 'Off';
end

seg_dir = ['./Results/SegmentedImages' sub_dir '/'];

seg_list = dir([seg_dir 'img_*.png']);
if(~isempty(seg_list))
    handles.ChooseDisplay.String{2} = 'Segmented';
    handles.ChooseDisplay.Value = 2;
    segIm = im2double(imread([seg_dir seg_list(1).name]));
    imshow(segIm);
end

output_dir = ['./Results/Tracking' sub_dir '/'];
if(exist([output_dir 'trajectory.mat']))
    handles.ChooseDisplay.String{3} = 'Tracking';
    handles.ChooseDisplay.Value = 3;
    handles.Progress.String = 'Tracking Complete.';
    set(handles.Text_Display, 'String', 'Tracking already completed for this folder');
    load([output_dir 'trajectory.mat']);
    load([output_dir 'norm_trajectory.mat']);
    load([output_dir 'CoM.mat']);
    
    handles.ViewerMode.Value = 1;
    handles.ViewTracking.Visible = 'On';
    handles.ViewTracking.Enable = 'On';
    
    startframe = find(trajectory, 1, 'first');
    endframe = length(CoM);
    
    % find body_length
    i = startframe;
    img = imread([seg_dir 'roi_' num2str(i) '.png']);
    img_norm = imtranslate(img, [255 - CoM(i,1), 255 - CoM(i,2)]);
    img_norm = imrotate(img_norm, CoM(i,3));
    img_norm = imcrop(img_norm, [size(img_norm,1)/2-150 size(img_norm,2)/2-150 300 300]);
    [Y,X] = find(img_norm);
    bodylength = max(Y(find(X==150)))-min(Y(find(X==150)));
    handles.BLText.Visible = 'On';
    handles.BodyLength.Visible = 'On';
    handles.BodyLength.String = num2str(round(bodylength));
    ResetImageSize_Callback(hObject, eventdata, handles);
    
    setappdata(hMainGui, 'start_frame',startframe);
    setappdata(hMainGui,'current_frame',startframe);
    setappdata(hMainGui, 'trajectory', trajectory);
    setappdata(hMainGui, 'norm_trajectory', norm_trajectory);
    setappdata(hMainGui, 'CoM', CoM);
    
    handles.VideoStart.Visible = 'On';
    handles.VideoEnd.Visible = 'On';
    handles.VideoStartFrame.Visible = 'On';
    handles.VideoEndFrame.Visible = 'On';
    handles.VideoStartFrame.String = num2str(startframe);
    handles.VideoEndFrame.String = num2str(endframe);

    nLegs = size(trajectory,2);
    if (nLegs > 6)
        colors = [1 0 0; 0 1 0; 0 0 1; 1 0 0.5; 1 1 0; 0 1 1; 1 0 1; 1 0.5 0];
        legs_id = {'L1', 'L2', 'L3','L4', 'R1', 'R2', 'R3','R4'};
    else
        colors = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 0 1 1; 1 0 1; 1 0 0.5; 1 0.5 0];
        legs_id = {'L1', 'L2', 'L3', 'R1', 'R2', 'R3'};
    end
    
    i = startframe;
    handles.text.String = ['Current frame: ' num2str(i)];
    handles.GoToDropdown.Value = i;
    
    handles.Tracking_Text.String = '';
    oriIm = im2double(imread([data_dir img_list(i).name]));
    imshow(oriIm);
    hold on;
    for kk = 1:size(trajectory,2)
        if(trajectory(i,kk,1)~=0)
            handles.Tracking_Text.String = [handles.Tracking_Text.String; sprintf('Leg %2s: %3d %3d', legs_id{kk}, trajectory(i,kk,1), trajectory(i,kk,2))]; 
            scatter(trajectory(i,kk,1),trajectory(i,kk,2),'w');
            text(trajectory(i,kk,1),trajectory(i,kk,2),legs_id{kk},'Color',[colors(kk,1) colors(kk,2) colors(kk,3)],'FontSize',14);
        else
            handles.Tracking_Text.String = [handles.Tracking_Text.String; sprintf('Leg %2s:        ', legs_id{kk})];
        end
    end
    hold off;
    
    drawnow
    
    % Update Display
    nLegs = 8 - handles.L1.Value - handles.L2.Value - handles.L3.Value - handles.R1.Value - handles.R2.Value - handles.R3.Value - handles.L4.Value - handles.R4.Value;
    if(length(str2num(handles.Tracking_Text.String(:,9:11))) < nLegs)
        set(handles.Missing, 'BackgroundColor', [1 0 0]); 
    else
        set(handles.Missing, 'BackgroundColor', [0 1 0]); 
    end
    
else
    setappdata(hMainGui, 'trajectory', []);
    setappdata(hMainGui, 'norm_trajectory', []);
    setappdata(hMainGui, 'CoM', []);
    handles.ViewerMode.Value = 0;
    handles.ViewTracking.Visible = 'Off';
    handles.ViewTracking.Enable = 'Off';
    handles.VideoStart.Visible = 'Off';
    handles.VideoEnd.Visible = 'Off';
    handles.VideoStartFrame.Visible = 'Off';
    handles.VideoEndFrame.Visible = 'Off';
end


%% --- Executes on button press in Segmentation.
function Segmentation_Callback(hObject, eventdata, handles)
% hObject    handle to Segmentation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hMainGui = getappdata(0, 'hMainGui');
data_dir = getappdata(hMainGui, 'data_dir');

set(handles.Text_Display, 'String', 'Starting segmentation.');
pause(1);
handles.axes1 = gcf;

pos_bs = strfind(data_dir,'Data');
sub_dir = data_dir(pos_bs(end)+length('Data'):length(data_dir));
data_dir = ['./Data' sub_dir '/'];
output_dir = ['./Results/SegmentedImages' sub_dir '/'];
if(~exist(output_dir))
    mkdir(output_dir);
end

if(~isempty(dir([data_dir '*.tif'])))
    fn_list = dir([data_dir '*.tif']);
else
    fn_list = dir([data_dir '*.bmp']);
end

I = imread([data_dir fn_list(1).name]);

if(~isempty(dir([data_dir 'Background/'])))
    if (exist([data_dir 'Background/Background.png']))
        ref_img = imread([data_dir 'Background/Background.png']);
        ref_img = double(ref_img);
        ref_img = padarray(ref_img,[20 20],'replicate');
        thres = 0.1;
        imshow(imcrop(ref_img,[21 21 size(I,2)-1 size(I,1)-1]),[]);
    else
        try
            background = dir([data_dir 'Background/']);
            ref_img = imread([data_dir 'Background/' background(3).name]);
            ref_img = double(ref_img);
            ref_img = padarray(ref_img,[20 20],'replicate');
            thres = 0.1;
            imshow(imcrop(ref_img,[21 21 size(I,2)-1 size(I,1)-1]),[]);
        catch
            [thres,~,~,ref_img,~] = video2background(data_dir, sub_dir);
            imshow(imcrop(ref_img,[21 21 size(I,2)-1 size(I,1)-1]),[]);
            if(~exist([data_dir 'Background']))
                mkdir([data_dir 'Background'])
            end
            imwrite(uint8(imcrop(ref_img,[21 21 size(I,2)-1 size(I,1)-1])),[data_dir 'Background/Background.png'], 'png');
        end
    end
else
    [thres,~,~,ref_img,~] = video2background(data_dir, sub_dir);
    imshow(imcrop(ref_img,[21 21 size(I,2)-1 size(I,1)-1]),[]);
    if(~exist([data_dir 'Background']))
        mkdir([data_dir 'Background'])
    end
    imwrite(uint8(imcrop(ref_img,[21 21 size(I,2)-1 size(I,1)-1])),[data_dir 'Background/Background.png'], 'png');
end
set(handles.Text_Display, 'String', 'This is the background averaged over the dataset.');
pause(5);

if(length(handles.ChooseDisplay.String) <2)
    handles.ChooseDisplay.String{2} = 'Segmented';
end

handles.ChooseDisplay.Value = 2;

params = setup_config_L('Drive');

%params = setup_lists(params);
%params = setup_directories_L(params);
params.border_skip_size = 20;
params.pos_samples_no = 30000;
params.neg_samples_no = 30000;
params.sample_size = [41 41]';

load_wl = handles.LoadClassifier.Value;
score_thres = str2num(handles.SegmentationConfidence.String);
if (load_wl)
    set(handles.Text_Display, 'String','Please select the trained classifier to be used.');
    [wl_file, wl_dir] = uigetfile(['./Results/Classifiers/','*.mat']);
    load([wl_dir wl_file]);
else
    if(~exist(['./Results/Classifiers' sub_dir '_Classifier.mat']))
        set(handles.Text_Display, 'String','Training the classifer over a random subset of the images.');
        sample_ratio = 20;
        tot_idx.pos = 0;
        tot_idx.neg = 0; 
        data = [];

        % Collect the positive and negative samples
        for i_img = 1 : floor(length(fn_list) / sample_ratio)

            I = imread([data_dir fn_list(i_img * sample_ratio).name]); 
            I = double(I);
            I = padarray(I,[20 20],'replicate');

            [pos_img,neg_img_body,neg_img_bkg] = leg_segment(I,ref_img,thres);
            show_img_output = repmat(I/255,[1 1 3]);
            show_img_output(:,:,1) = show_img_output(:,:,1) + pos_img;
            show_img_output(:,:,3) = show_img_output(:,:,3) + neg_img_body;
            imshow(show_img_output,[]);

            data.train.imgs.X{i_img,1} = I / 255;
            data.train.imgs.X{i_img,2} = (ref_img - I) / 255;

            gt_img = zeros(size(I));  
            gt_img(pos_img > 0) = 1;    
            gt_img(neg_img_body > 0) = -1;      
            gt_img(neg_img_bkg > 0) = -1;
            data.train.gts.X{i_img} = gt_img;

            % ------------------------------

            border_img = zeros(size(I));
            border_img(1:params.border_skip_size,:) = 1;
            border_img(:,1:params.border_skip_size) = 1;
            border_img((end-params.border_skip_size):end,:) = 1;
            border_img(:,(end-params.border_skip_size):end) = 1;

            pos_sampling = find((pos_img) & (~border_img) & (~neg_img_body)); 
            neg_body_sampling = find(neg_img_body & (~border_img) & (~pos_img));
            neg_bkg_sampling = find(neg_img_bkg & (~border_img));

            % sample counts
            imgs_no = floor(length(fn_list) / sample_ratio);
            npos_img = ceil(params.pos_samples_no / imgs_no);
            nneg_img = ceil(params.neg_samples_no / imgs_no);

            % getting a random subset?
            neg_sampling = neg_body_sampling(randi(length(neg_body_sampling),[nneg_img,1]));
            neg_sampling = [neg_sampling; neg_bkg_sampling(randi(length(neg_bkg_sampling),[nneg_img,1]))];

            npos = length(pos_sampling);        
            nneg = length(neg_sampling);

            good_sampling_idx{i_img}.pos = pos_sampling(randi(npos,[ceil(npos_img),1]));        
            good_sampling_idx{i_img}.neg = neg_sampling(randi(nneg,[ceil(nneg_img),1]));
            good_sampling_idx{i_img}.pos_val = pos_img(good_sampling_idx{i_img}.pos);        
            good_sampling_idx{i_img}.neg_val = pos_img(good_sampling_idx{i_img}.neg);

            tot_idx.pos = tot_idx.pos+length(good_sampling_idx{i_img}.pos);        
            tot_idx.neg = tot_idx.neg+length(good_sampling_idx{i_img}.neg);
        end

        data.train.imgs.idxs = 1:2;
        data.train.imgs.sub_ch_no = 2;
        samples_idx = [];
        
        % Preparing the samples
        for i_img = 1 : floor(length(fn_list) / sample_ratio)

            I = imread([data_dir fn_list(i_img * sample_ratio).name]);
            I = padarray(I,[20 20],'replicate');

            samp_idx = [good_sampling_idx{i_img}.pos ; good_sampling_idx{i_img}.neg];

            samp_idx_2D = zeros(length(samp_idx),2);
            [samp_idx_2D(:,1),samp_idx_2D(:,2)] = ind2sub(size(I),samp_idx);

            labels = [good_sampling_idx{i_img}.pos_val ; good_sampling_idx{i_img}.neg_val];
            labels = double(labels);
            labels(labels < 0.5) = -1;

            samples_idx_img = zeros(length(samp_idx),4);
            samples_idx_img(:,1) = i_img;
            samples_idx_img(:,2:3) = samp_idx_2D;
            samples_idx_img(:,4) = labels;

            samples_idx = cat(1,samples_idx,samples_idx_img);

        end

        samples_idx = samples_idx(1 : (params.pos_samples_no + params.neg_samples_no),:);

        params.pos_samples_no = sum(samples_idx(:,end)==1);   
        params.neg_samples_no = sum(samples_idx(:,end)==-1);  

        params.T1_size = round((params.pos_samples_no+params.neg_samples_no)/3);
        params.T2_size = round((params.pos_samples_no+params.neg_samples_no)/3);

        params.pos_to_sample_no = params.T1_size;
        params.neg_to_sample_no = params.T1_size;
        params.wl_no = 100;

        % Train the classifier here
        weak_learners = train_boost_general(params,data,samples_idx);
        sub_dir_pos = strfind(sub_dir,'/');
        if(~exist(['./Results/Classifiers/' sub_dir(1:sub_dir_pos(end))]))
            mkdir(['./Results/Classifiers/' sub_dir(1:sub_dir_pos(end))]);
        end
        save (['./Results/Classifiers/' sub_dir '_Classifier.mat'],'weak_learners');
    else
        load(['./Results/Classifiers/' sub_dir '_Classifier.mat']);
    end    
end

set(handles.Text_Display, 'String','Training completed, now applying the classifier to all images.');
fprintf('Output directory: %s\n', output_dir);

% Applying the classifier
sample_ratio = 1;
sec_no = 10;

for i_sec = 1 : ceil(length(fn_list) / sample_ratio / sec_no)
    X = [];
    if (i_sec < ceil(length(fn_list) / sample_ratio / sec_no))
        imgs_sec = 1 + (i_sec - 1) * sec_no : sample_ratio :(i_sec) * sec_no;
    else
        imgs_sec = 1 + (i_sec - 1) * sec_no : sample_ratio : length(fn_list);
    end
      
    for i_img = 1 : length(imgs_sec)
        I = imread([data_dir fn_list(imgs_sec(i_img)).name]);
        I = double(I);
        I = padarray(I,[20 20],'replicate');
        
        roi_img = (max(ref_img - I(:,:,1),0) ./ ref_img) > thres;
        roi_img = bwareaopen(roi_img, 200);
        
        %neg_img{i_img} = leg_segment_clean(I,roi_img);
        
        bkg_sub = (ref_img - I);        
        I(:,:,2) = bkg_sub;

        I = I / 255;    
        
        roi_images{i_img} = roi_img;
        
        X{i_img,1} = I(:,:,1);
        X{i_img,2} = I(:,:,2);
    end
    
    score_images = batch_evaluate_boost_images(X,params,weak_learners,roi_images);
    
    for i_img = 1 : length(imgs_sec)
        I = imread([data_dir fn_list(imgs_sec(i_img)).name]);        
        I = double(I)/255;
        I = padarray(I,[20 20],'replicate');       
        
        score_img = score_images{i_img};      
        est_leg = score_img > score_thres;
        %est_leg = bwdist(est_leg)<3 .*roi_images{i_img};
        %est_leg = est_leg & ~ neg_img{i_img} & roi_images{i_img};
        est_leg = filter_small_comp(est_leg,5);
    
        show_img_output = repmat(I,[1 1 3]);
        show_img_output(:,:,1) = show_img_output(:,:,1) + est_leg;
        
        imwrite(imcrop(show_img_output,[21 21 size(I,2)-41 size(I,1)-41]),[output_dir 'img_' num2str(imgs_sec(i_img)) '.png'],'png');
        imwrite(imcrop(roi_images{i_img},[21 21 size(I,2)-41 size(I,1)-41]),[output_dir 'roi_' num2str(imgs_sec(i_img)) '.png'],'png');
        if (i_img == 1)
            imshow(imcrop(show_img_output,[21 21 size(I,2)-41 size(I,1)-41]));
            handles.GoToDropdown.Value = (i_sec-1)*10 + i_img;
            handles.Progress.String = ['Segmentation Progress: ' num2str(min(i_sec*10,length(fn_list))) '/' num2str(length(fn_list))];
            pause(0.5);
        end
    end
end

handles.Text_Display.String = 'Segmentation complete, proceeding with tracking.';
pause(5);
Tracking_Callback(hObject, eventdata, handles);


%% --- Executes on button press in Tracking.
function Tracking_Callback(hObject, eventdata, handles)
% hObject    handle to Tracking (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hMainGui = getappdata(0, 'hMainGui');
data_dir = getappdata(hMainGui, 'data_dir');
if(isempty(getappdata(hMainGui, 'trajectory')))
    setappdata(hMainGui, 'CoM', []);
    setappdata(hMainGui, 'trajectory', []);
    setappdata(hMainGui, 'norm_trajectory', []);
    setappdata(hMainGui, 'start_frame', []);
end
CoM = getappdata(hMainGui, 'CoM');
trajectory = getappdata(hMainGui, 'trajectory');
norm_trajectory = getappdata(hMainGui, 'norm_trajectory');

pos_bs = strfind(data_dir,'Data');
sub_dir = data_dir(pos_bs(end)+length('Data'):length(data_dir));
seg_dir = ['./Results/SegmentedImages' sub_dir '/'];
data_dir = ['./Data' sub_dir '/'];
if(~isempty(dir([data_dir '*.tif'])))
    img_list = dir([data_dir '*.tif']);
else
    img_list = dir([data_dir '*.bmp']);
end
output_dir = ['./Results/Tracking/' sub_dir '/'];
if(~exist(output_dir))
    mkdir(output_dir);
end

if(length(handles.ChooseDisplay.String) <3)
    handles.ChooseDisplay.String{3} = 'Tracking';
end

removed_legs = handles.L1.Value+handles.L2.Value+handles.L3.Value+handles.L4.Value+handles.R1.Value+handles.R2.Value+handles.R3.Value+handles.R4.Value;

nLegs = 8 - handles.L4.Value - handles.R4.Value;
if (nLegs > 6)
    colors = [1 0 0; 0 1 0; 0 0 1; 1 0 0.5; 1 1 0; 0 1 1; 1 0 1; 1 0.5 0];
    legs_id = {'L1', 'L2', 'L3','L4', 'R1', 'R2', 'R3','R4'};
else
    colors = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 0 1 1; 1 0 1; 1 0 0.5; 1 0.5 0];
    legs_id = {'L1', 'L2', 'L3', 'R1', 'R2', 'R3'};
end

se = strel('disk',4);
max_distance = 20;

if (isempty(seg_dir))
    set(handles.Text_Display, 'String', 'Please ensure that segmentation has been completed first.');
else
    set(handles.Text_Display, 'String', 'Performing legs tracking.');
    pause(1);
    handles.axes1 = gcf;
    
    drawnow
    
    end_frame = length(dir([seg_dir 'img_*.png']));
    handles.ChooseDisplay.Value = 3;
    
    %Manually initiate tracking by labeling leg tips on chosen frame 
    if(get(handles.ManualInitiate,'Value'))
        setappdata(hMainGui, 'CoM', []);
        setappdata(hMainGui, 'trajectory', []);
        setappdata(hMainGui, 'norm_trajectory', []);
        
        i = handles.GoToDropdown.Value;
        Imroi = imread([seg_dir 'roi_' num2str(i) '.png']);
        handles.Tracking_Text.String = '';
        setappdata(hMainGui, 'current_frame', i);
        setappdata(hMainGui, 'start_frame', i);
        handles.GoToDropdown.String = i:i;
        
        Imroi = imerode(Imroi,se);
        [locY,locX] = find(Imroi > 0);

        %Centre_of_Mass
        [centre_of_mass,theta] = findCoM(locX,locY);

        CoM(i, :) = [centre_of_mass, theta];
        setappdata(hMainGui, 'CoM' , CoM);
        oriIm = im2double(imread([data_dir img_list(i).name]));
        imshow(oriIm);
        handles.Tracking.String = 'Initial';
        
    %Automatically initiate tracking
    elseif(isempty(getappdata(hMainGui, 'start_frame')));
        for i = 1: end_frame
            clear connComp;
            % Load necessary images
            Im = imread([seg_dir 'img_' num2str(i) '.png']);
            Imroi = imread([seg_dir 'roi_' num2str(i) '.png']);

            % Fix image/segmentation irregularities
            Imwork = double(Im(:,:,1) == 255) .* Imroi;
            connComp = bwconncomp(Imwork);
            
            %Find a frame with correct number of legs (each leg component
            %must contain at least 35 pixels)
            if (connComp.NumObjects == 8 - removed_legs && min(cellfun(@length,connComp.PixelIdxList)) > 35)
                handles.text.String = ['Current frame: ' num2str(i)];

                Imroi = imerode(Imroi,se);
                [locY,locX] = find(Imroi > 0);

                %Centre_of_Mass
                [centre_of_mass,theta] = findCoM(locX,locY);

                locX_norm = locX - centre_of_mass(1);
                locY_norm = locY - centre_of_mass(2);
                points = [cosd(theta) sind(theta); -sind(theta) cosd(theta)]* [locX_norm'; locY_norm'];
                
                %Attempt to ensure that the head is in the positive
                %y-direction by comparing which side has more number of
                %pixel points
                if (length(find(points(2,:) > 0)) < length(find(points(2,:) < 0)))
                    theta = mod(theta + 180,360);
                end
                
                for cc = 1 : 8 - removed_legs
                    legpixels = connComp.PixelIdxList{cc};

                    imleg = zeros(size(Imwork));
                    imleg(legpixels) = 1;

                    [skr,~] = skeleton(imleg);
                    skr_start = 3;
                    [~,exy,~] = anaskel(bwmorph(skr > skr_start,'skel',Inf));
                    
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
                if(length(leftlegs_idx)~=4 - handles.L1.Value - handles.L2.Value - handles.L3.Value - handles.L4.Value)
                    continue;
                end
                rightlegs_idx = find(raw_normtips(:,1)>0);
                if(length(rightlegs_idx)~=4 - handles.R1.Value - handles.R2.Value - handles.R3.Value - handles.R4.Value)
                    continue;
                end

                [~,left_sort] = sort(raw_normtips(leftlegs_idx,2));

                leg_counter = 0;
                for kk = 1 : nLegs/2
                    if (eval(['handles.L' num2str(kk) '.Value']))
                        leg_counter = leg_counter+1;
                        continue;
                    end
                    trajectory(i, kk, :) = raw_tips(leftlegs_idx(left_sort(kk-leg_counter)),:);
                    norm_trajectory(i, kk, :) = raw_normtips(leftlegs_idx(left_sort(kk-leg_counter)),:);
                end

                [~,right_sort] = sort(raw_normtips(rightlegs_idx,2));
                
                leg_counter = 0;
                for kk = 1 :nLegs/2
                    if (eval(['handles.R' num2str(kk) '.Value']))
                        leg_counter = leg_counter+1;
                        continue;
                    end
                    trajectory(i, kk+nLegs/2, :) = raw_tips(rightlegs_idx(right_sort(kk-leg_counter)),:);
                    norm_trajectory(i, kk+nLegs/2, :) = raw_normtips(rightlegs_idx(right_sort(kk-leg_counter)),:);
                end
                
                CoM(i, :) = [centre_of_mass, theta];

                handles.Tracking_Text.String = '';

                setappdata(hMainGui, 'start_frame', i);
                setappdata(hMainGui, 'current_frame', i);

                oriIm = im2double(imread([data_dir img_list(i).name]));
                imshow(oriIm);
                hold on;
                
                % find body_length
                img = imread([seg_dir 'roi_' num2str(i) '.png']);
                img_norm = imtranslate(img, [255 - CoM(i,1), 255 - CoM(i,2)]);
                img_norm = imrotate(img_norm, CoM(i,3));
                img_norm = imcrop(img_norm, [size(img_norm,1)/2-150 size(img_norm,2)/2-150 300 300]);
                [Y,X] = find(img_norm);
                bodylength = max(Y(find(X==150)))-min(Y(find(X==150)));
                handles.BLText.Visible = 'On';
                handles.BodyLength.Visible = 'On';
                handles.BodyLength.String = num2str(round(bodylength));
                ResetImageSize_Callback(hObject, eventdata, handles);
                
                for kk = 1:nLegs
                    if (eval(['handles.' legs_id{kk} '.Value']))
                        handles.Tracking_Text.String = [handles.Tracking_Text.String; sprintf('Leg %2s:        ', legs_id{kk})];
                    end
                    scatter(trajectory(i,kk,1),trajectory(i,kk,2),'w');
                    text(trajectory(i,kk,1),trajectory(i,kk,2),legs_id{kk},'Color',[colors(kk,1) colors(kk,2) colors(kk,3)],'FontSize',14);
                    handles.Tracking_Text.String = [handles.Tracking_Text.String; sprintf('Leg %2s: %3d %3d', legs_id{kk}, trajectory(i,kk,1), trajectory(i,kk,2))];
                end
                hold off;
                drawnow
                
                for k = i-1 : -1 : 1
                    Im = (imread([seg_dir 'img_' num2str(k) '.png']));
                    Imroi = imread([seg_dir 'roi_' num2str(k) '.png']);

                    % Fix image/segmentation irregularities
                    Imwork = double(Im(:,:,1) == 255) .* Imroi;
                    if (nLegs>6)
                        Imwork = Imwork - (bwdist(imerode(Imroi,strel('disk',6)))<15);
						Imwork = bwareaopen(Imwork>0,5);
                    else
                        Imwork = Imwork - (bwdist(imerode(Imroi,se))<10);
						Imwork = bwareaopen(Imwork>0,3);
                    end
                    
                    [locY,locX] = find(imerode(Imroi,se) > 0);
                    
                    [centre_of_mass,theta] = findCoM(locX,locY);

                    if (abs(theta - CoM(k+1,3)) > 90 && abs(theta - CoM(k+1,3)) < 270)
                        theta = mod(theta + 180,360);
                    end

                    CoM(k, :) = [centre_of_mass, theta];

                    clear raw_tips;
                    clear raw_normtips;
                    
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
                        if(eval(['handles.' legs_id{unassigned_idx(kk)} '.Value']))
                            continue;
                        end
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
    end    
    
    if (strcmp(handles.Tracking.String, 'Tracking') || strcmp(handles.Tracking.String, 'Resume'))
        handles.Tracking.String = 'Pause';
    else
        if (~strcmp(handles.Tracking.String, 'Initial'))
            handles.Tracking.String = 'Resume';
        end
    end
    
    skip = 1;
   
    while (strcmp(handles.Tracking.String, 'Pause'))
        setappdata(hMainGui, 'current_frame', getappdata(hMainGui, 'current_frame') + skip);
        i = getappdata(hMainGui, 'current_frame');
        if (i>end_frame)
            setappdata(hMainGui, 'current_frame', end_frame);
            break;
        end
        
        handles.Progress.String = ['Tracking Progress: ' num2str(i) '/' num2str(length(img_list))];
        handles.text.String = ['Current frame: ' num2str(i)];

        % Load necessary images
        Im = (imread([seg_dir 'img_' num2str(i) '.png']));
        Imroi = imread([seg_dir 'roi_' num2str(i) '.png']);
    
        % Fix image/segmentation irregularities
        Imwork = double(Im(:,:,1) == 255) .* Imroi;
        [locY,locX] = find(imerode(Imroi,se) > 0);
        
        if (nLegs>6)
            Imwork = Imwork - (bwdist(imerode(Imroi,strel('disk',6)))<15);
						Imwork = bwareaopen(Imwork>0,5);
        else
            Imwork = Imwork - (bwdist(imerode(Imroi,se))<10);
			Imwork = bwareaopen(Imwork>0,3);
        end
        
        %Centre_of_Mass
        [centre_of_mass,theta] = findCoM(locX,locY);
        
        if (abs(theta - CoM(i-skip,3)) > 90 && abs(theta - CoM(i-skip,3)) < 270)
            theta = mod(theta + 180,360);
        end

        CoM(i, :) = [centre_of_mass, theta];

        clear raw_tips;
        clear raw_normtips;
        
        [raw_tips, raw_normtips] = findtips(Imwork, Imroi, locX, locY, centre_of_mass, theta);
        
        x_t1 = raw_normtips;
        x_t0 = reshape(norm_trajectory(i-1,:),[nLegs 2]);
        target_indices = hungarianlinker(x_t0, x_t1, max_distance);

        for kk = 1:length(target_indices)
            if (target_indices(kk) > 0)
                norm_trajectory(i, kk, :) = raw_normtips(target_indices(kk),:);
                trajectory(i, kk, :) = raw_tips(target_indices(kk),:);
            else
                norm_trajectory(i, kk, :) = 0;
                trajectory(i, kk, :) = 0;
            end
        end

        leftoveridx = find(hungarianlinker(raw_normtips, x_t0, max_distance)==-1);
        x_t1 = raw_normtips(leftoveridx,:);
        x_t2 = raw_tips(leftoveridx,:);

        last_seen_tips = reshape(norm_trajectory(i,:),[nLegs 2]);
        unassigned_idx = find(last_seen_tips(:,1) == 0);
        for kk = 1:length(unassigned_idx)
            if(eval(['handles.' legs_id{unassigned_idx(kk)} '.Value']))
                continue;
            end
            idx = find(norm_trajectory(1:i,unassigned_idx(kk),1),1,'last');
            last_seen_tips(unassigned_idx(kk),:) = reshape(norm_trajectory(idx,unassigned_idx(kk),:), [1 2]);
        end

        if (~isempty(unassigned_idx) && ~isempty(x_t1))
                C = hungarianlinker(last_seen_tips(unassigned_idx,:), x_t1, 1.25*max_distance);
                for l = 1:length(C)
                    if (C(l) ~= -1)
                        norm_trajectory(i,unassigned_idx(l), :) = x_t1(C(l), :);
                        trajectory(i,unassigned_idx(l), :) = x_t2(C(l), :);
                    end
                end
        end
        
        handles.GoToDropdown.Value = i;
        handles.Tracking_Text.String = '';
        oriIm = im2double(imread([data_dir img_list(i).name]));
        imshow(oriIm);
        hold on;
        for kk = 1:nLegs
            if (eval(['handles.' legs_id{kk} '.Value']))
                trajectory(i,kk,1) = 0;
                trajectory(i,kk,2) = 0;
                norm_trajectory(i,kk,1) = 0;
                norm_trajectory(i,kk,2) = 0;
                handles.Tracking_Text.String = [handles.Tracking_Text.String; sprintf('Leg %2s:        ', legs_id{kk})];
                continue;
            end
            if(trajectory(i,kk,1)~=0)
                handles.Tracking_Text.String = [handles.Tracking_Text.String; sprintf('Leg %2s: %3d %3d', legs_id{kk}, trajectory(i,kk,1), trajectory(i,kk,2))]; 
                scatter(trajectory(i,kk,1),trajectory(i,kk,2),'w');
                text(trajectory(i,kk,1),trajectory(i,kk,2),legs_id{kk},'Color',[colors(kk,1) colors(kk,2) colors(kk,3)],'FontSize',14);
            else
                handles.Tracking_Text.String = [handles.Tracking_Text.String; sprintf('Leg %2s:        ', legs_id{kk})];
            end
        end
        hold off;
        drawnow
        
        % Update Display
        nLegs = 8 - handles.L1.Value - handles.L2.Value - handles.L3.Value - handles.R1.Value - handles.R2.Value - handles.R3.Value - handles.L4.Value - handles.R4.Value;
        if(length(str2num(handles.Tracking_Text.String(:,9:11))) < nLegs)
            set(handles.Missing, 'BackgroundColor', [1 0 0]); 
        else
            set(handles.Missing, 'BackgroundColor', [0 1 0]); 
        end
        
        setappdata(hMainGui, 'CoM' , CoM);
        setappdata(hMainGui, 'trajectory' , trajectory);
        setappdata(hMainGui, 'norm_trajectory' , norm_trajectory);
        
        handles.GoToDropdown.String = 1:find(CoM(:,1),1,'last');
    end
    
    if (getappdata(hMainGui, 'current_frame') == end_frame)
        save([output_dir 'CoM.mat'],'CoM');
        save([output_dir 'trajectory.mat'],'trajectory');
        save([output_dir 'norm_trajectory.mat'],'norm_trajectory');
        csvwrite([output_dir 'CoM.csv'],CoM);
        csvwrite([output_dir 'trajectory.csv'],trajectory);
        csvwrite([output_dir 'norm_trajectory.csv'],norm_trajectory);
        handles.Tracking.String = 'Complete';
        handles.Text_Display.String = 'Tracking Complete';
        handles.VideoStart.Visible = 'On';
        handles.VideoEnd.Visible = 'On';
        handles.VideoStartFrame.Visible = 'On';
        handles.VideoEndFrame.Visible = 'On';
        handles.VideoStartFrame.String = num2str(getappdata(hMainGui, 'start_frame'));
        handles.VideoEndFrame.String = num2str(end_frame);
    end
end


%% --- Mouse cursor location in axis
function fh_wbmfcn(varargin)
% WindowButtonMotionFcn for the figure.
S = varargin{3};  % Get the structure.
F = get(S.fh,'currentpoint');  % The current point w.r.t the figure.
% Figure out of the current point is over the axes or not -> logicals.
tf1 = (S.fh.Position(3)-512)/2 <= F(1) && F(1) <= (S.fh.Position(3)+512)/2;
tf2 = (S.fh.Position(4)-512) <= F(2) && F(2) <= S.fh.Position(4);

if tf1 && tf2
    % Calculate the current point w.r.t. the axes.
    Cx = round(F(1) - (S.fh.Position(3)-512)/2);
    Cy = round(S.fh.Position(4) - F(2));
    % The mouse location is displayed at handles.Position_Text.
    set(S.tx,'str',num2str([Cx,Cy],'%3.0f   %3.0f'))
end


%% --- Double click to label tip position of a leg
function clicker(h,~,handles)
    hMainGui = getappdata(0, 'hMainGui');
    nLegs = 8  - handles.L4.Value - handles.R4.Value;
    if (nLegs > 6)
        colors = [1 0 0; 0 1 0; 0 0 1; 1 0 0.5; 1 1 0; 0 1 1; 1 0 1; 1 0.5 0];
        legs_id = {'L1', 'L2', 'L3','L4', 'R1', 'R2', 'R3','R4'};
        switch get(get(handles.Labelling_Panel, 'SelectedObject'), 'Tag');
            case 'Leg_L1', id = 1;
            case 'Leg_L2', id = 2;
            case 'Leg_L3', id = 3;
            case 'Leg_L4', id = 4;
            case 'Leg_R1', id = 5;
            case 'Leg_R2', id = 6;
            case 'Leg_R3', id = 7;
            case 'Leg_R4', id = 8;
        end
    else
        colors = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 0 1 1; 1 0 1; 1 0 0.5; 1 0.5 0];
        legs_id = {'L1', 'L2', 'L3', 'R1', 'R2', 'R3'};
        switch get(get(handles.Labelling_Panel, 'SelectedObject'), 'Tag');
            case 'Leg_L1', id = 1;
            case 'Leg_L2', id = 2;
            case 'Leg_L3', id = 3;
            case 'Leg_R1', id = 4;
            case 'Leg_R2', id = 5;
            case 'Leg_R3', id = 6;
        end
    end
    
    F = get(h, 'currentpoint');
    if (strcmp(get(h, 'selectiontype'),'open') && (handles.figMain.Position(3)-512)/2 <= F(1) && F(1) <= (handles.figMain.Position(3)+512)/2 && (handles.figMain.Position(4)-512) <= F(2) && F(2) <= handles.figMain.Position(4))
        Cx = round(F(1) - (handles.figMain.Position(3)-512)/2);
        Cy = round(handles.figMain.Position(4) - F(2));
        if (strcmp(handles.Annotation.String,'Save and Exit'))
            handles.Tracking_Text.String(id,:) = sprintf('Leg %2s: %3d %3d', legs_id{id}, Cx, Cy);
            for i = 1: length(handles.figMain.CurrentAxes.Children)-2
               if(strcmp(handles.figMain.CurrentAxes.Children(i).Type,'text') && strcmp(handles.figMain.CurrentAxes.Children(i).String,legs_id{id}))
                   delete(handles.figMain.CurrentAxes.Children(i:i+1));
               end
            end
            hold on
            scatter(Cx,Cy,'w');
            text(Cx,Cy,legs_id{id},'Color',[colors(id,1) colors(id,2) colors(id,3)],'FontSize',14);
            hold off
        end
    end

    
function update_display(hObject, eventdata, handles)
nLegs = 8 - handles.L1.Value - handles.L2.Value - handles.L3.Value - handles.R1.Value - handles.R2.Value - handles.R3.Value - handles.L4.Value - handles.R4.Value;
if(length(str2num(handles.Tracking_Text.String(:,9:11))) < nLegs)
    set(handles.Missing, 'BackgroundColor', [1 0 0]); 
else
    set(handles.Missing, 'BackgroundColor', [0 1 0]); 
end


%% --- Executes when selected object is changed in Labelling_Panel.
function Labelling_Panel_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in Labelling_Panel 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (strcmp(handles.Tracking.String,'Resume'))
switch get(get(handles.Labelling_Panel, 'SelectedObject'), 'Tag');
    case 'Leg_L1', set(handles.Text_Display,'String','Selecting leg L1');
    case 'Leg_L2', set(handles.Text_Display,'String','Selecting leg L2');
    case 'Leg_L3', set(handles.Text_Display,'String','Selecting leg L3');
    case 'Leg_R1', set(handles.Text_Display,'String','Selecting leg R1');
    case 'Leg_R2', set(handles.Text_Display,'String','Selecting leg R2');
    case 'Leg_R3', set(handles.Text_Display,'String','Selecting leg R3');
    case 'Leg_L4', set(handles.Text_Display,'String','Selecting leg L4');
    case 'Leg_R4', set(handles.Text_Display,'String','Selecting leg R4');
end
end


%% --- Executes on button press in Video.
function Video_Callback(hObject, eventdata, handles)
% hObject    handle to Video (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hMainGui = getappdata(0, 'hMainGui');
data_dir = getappdata(hMainGui, 'data_dir');
pos_bs = strfind(data_dir,'Data');
sub_dir = data_dir(pos_bs(end)+length('Data'):length(data_dir));
tracking_dir = ['./Results/Tracking' sub_dir '/'];
if (exist(['./Results/Tracking' sub_dir '/trajectory.mat']))
    figure
    Video_base(data_dir,30,str2num(handles.Video_Step.String),str2num(handles.VideoStartFrame.String),str2num(handles.VideoEndFrame.String));
else
    handles.Text_Display.String = 'Please ensure that tracking has been completed first.';
end


%% --- Executes on selection change in GoToDropdown.
function GoToDropdown_Callback(hObject, eventdata, handles)
% hObject    handle to GoToDropdown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns GoToDropdown contents as cell array
%        contents{get(hObject,'Value')} returns selected item from GoToDropdown
hMainGui = getappdata(0, 'hMainGui');
i = handles.GoToDropdown.Value;
setappdata(hMainGui,'current_frame',i);
handles.text.String = ['Current frame: ' num2str(i)];
data_dir = getappdata(hMainGui, 'data_dir');
pos_bs = strfind(data_dir,'Data');
sub_dir = data_dir(pos_bs(end)+length('Data'):length(data_dir));
data_dir = ['./Data' sub_dir '/'];
if(~isempty(dir([data_dir '*.tif'])))
    img_list = dir([data_dir '*.tif']);
else
    img_list = dir([data_dir '*.bmp']);
end
oriIm = im2double(imread([data_dir img_list(i).name]));
imshow(oriIm);

display = handles.ChooseDisplay.Value;

switch display
    case 1
        handles.GoToDropdown.String = 1:length(img_list);
        imshow(oriIm);
    case 2
        seg_dir = ['./Results/SegmentedImages' sub_dir '/'];
        seg_list = dir([seg_dir 'img_*.png']);
        handles.GoToDropdown.String = 1:length(seg_list);
        imshow([seg_dir 'img_' num2str(i) '.png']);
    case 3
        trajectory = getappdata(hMainGui, 'trajectory');
        startframe = 1;
        handles.Tracking_Text.String = '';
        
        i = handles.GoToDropdown.Value + startframe - 1;
        setappdata(hMainGui,'current_frame',i);
        handles.text.String = ['Current frame: ' num2str(i)];
        
        oriIm = im2double(imread([data_dir img_list(i).name]));
        imshow(oriIm);
        
        handles.GoToDropdown.String = 1:length(trajectory);
        
        nLegs = size(trajectory,2);
        if (nLegs > 6)
            colors = [1 0 0; 0 1 0; 0 0 1; 1 0 0.5; 1 1 0; 0 1 1; 1 0 1; 1 0.5 0];
            legs_id = {'L1', 'L2', 'L3','L4', 'R1', 'R2', 'R3','R4'};
        else
            colors = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 0 1 1; 1 0 1; 1 0 0.5; 1 0.5 0];
            legs_id = {'L1', 'L2', 'L3', 'R1', 'R2', 'R3'};
        end
        
        hold on;
        for kk = 1:nLegs
            if(eval(['handles.' legs_id{kk} '.Value']))
                handles.Tracking_Text.String = [handles.Tracking_Text.String; sprintf('Leg %2s:        ', legs_id{kk})];
                continue;
            end
            if(trajectory(i,kk,1)~=0)
                handles.Tracking_Text.String = [handles.Tracking_Text.String; sprintf('Leg %2s: %3d %3d', legs_id{kk}, trajectory(i,kk,1), trajectory(i,kk,2))]; 
                scatter(trajectory(i,kk,1),trajectory(i,kk,2),'w');
                text(trajectory(i,kk,1),trajectory(i,kk,2),legs_id{kk},'Color',[colors(kk,1) colors(kk,2) colors(kk,3)],'FontSize',14);
            else
                handles.Tracking_Text.String = [handles.Tracking_Text.String; sprintf('Leg %2s:        ', legs_id{kk})];
            end
        end
        hold off;
        drawnow
        
        % Update Display
        nLegs = 8 - handles.L1.Value - handles.L2.Value - handles.L3.Value - handles.R1.Value - handles.R2.Value - handles.R3.Value - handles.L4.Value - handles.R4.Value;
        if(length(str2num(handles.Tracking_Text.String(:,9:11))) < nLegs)
            set(handles.Missing, 'BackgroundColor', [1 0 0]); 
        else
            set(handles.Missing, 'BackgroundColor', [0 1 0]); 
        end
end


%% --- Executes on button press in BatchProcess.
function BatchProcess_Callback(hObject, eventdata, handles)
% hObject    handle to BatchProcess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hMainGui = getappdata(0, 'hMainGui');
prompt = {'Input 0 for segmentation only or 1 to include tracking'};
dlg_title = 'Options for batch process';
num_lines = 1;
defaultans = {'0'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
if(isempty(answer))
    return;
end
answer = ~(answer{1}=='0');
data_dir = uipickfiles;
try
    if(data_dir==0)
        return;
    end
catch
    
end
if(isempty(data_dir))
    return;
end
for i = 1 : length(data_dir)
    ProcessFolder(data_dir{i}, answer, handles);
end


function ProcessFolder(data_dir, answer, handles)
score_thres = str2num(handles.SegmentationConfidence.String);
pos_bs = strfind(data_dir,'Data');
sub_dir = data_dir(pos_bs(end)+length('Data'):length(data_dir));
if(exist(data_dir) && isempty(dir([data_dir '/*.tif'])) && isempty(dir([data_dir '/*.bmp'])))
    dir_fn = dir(data_dir);
    for i_dir = 1:length(dir_fn)-2
        if (dir_fn(i_dir+2).isdir)
            ProcessFolder([data_dir '/' dir_fn(i_dir+2).name], answer, handles);
            pause(1);
        end
    end
else
    data_dir = ['./Data' sub_dir];
    
    handles.Text_Display.String = ['Performing segmentation for the folder: ' data_dir];
    fprintf('Performing segmentation for the folder: %s.\n', data_dir);
    pause(1);
    try
        Segmentation(data_dir,score_thres,0);
    catch MExc
        fprintf('An error occured during segmentation for the folder: %s.\n', data_dir);
        disp('Execution will continue.');
        disp(MExc.message);
    end
    handles.Text_Display.String = ['Completed segmentation for the folder: ' data_dir];
    fprintf('Completed segmentation for the folder: %s.\n', data_dir);
    pause(1)
    if (answer)
        handles.Text_Display.String = ['Performing tracking for the folder: ' data_dir];
        fprintf('Performing tracking for the folder: %s.\n', data_dir);
        pause(1);
        try    
            Tracking_Base(data_dir);
        catch MExc
            fprintf('An error occured during tracking for the folder: %s.\n', data_dir);
            disp('Execution will continue.');
            disp(MExc.message);
        end
        handles.Text_Display.String = ['Completed tracking for the folder: ' data_dir];
        fprintf('Completed tracking for the folder: %s.\n', data_dir);
        pause(1);
    end
end


%% --- Change display mode
function ChooseDisplay_Callback(hObject, eventdata, handles)
% hObject    handle to ChooseDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ChooseDisplay contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ChooseDisplay
GoToDropdown_Callback(hObject, eventdata, handles);


%% --- Change start frame of make video
function VideoStartFrame_Callback(hObject, eventdata, handles)
% hObject    handle to VideoStartFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hMainGui = getappdata(0, 'hMainGui');
if (str2num(hObject.String) < getappdata(hMainGui,'start_frame'))
    hObject.String = num2str(getappdata(hMainGui,'start_frame'));
end


%% --- Change end frame of make video
function VideoEndFrame_Callback(hObject, eventdata, handles)
% hObject    handle to VideoEndFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hMainGui = getappdata(0, 'hMainGui');
CoM = getappdata(hMainGui, 'CoM');
if (str2num(hObject.String) > length(CoM))
    hObject.String = num2str(length(CoM));
end
if (str2num(hObject.String) < str2num(handles.VideoStartFrame.String))
    hObject.String = num2str(str2num(handles.VideoStartFrame.String)+str2num(handles.Video_Step.String));
end


%% --- Change segmentation threshold via slider
function SlideThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to SlideThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.SegmentationConfidence.String = num2str(handles.SlideThreshold.Value);


%% --- Change segmentation threshold via numerical input
function SegmentationConfidence_Callback(hObject, eventdata, handles)
% hObject    handle to SegmentationConfidence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SegmentationConfidence as text
%        str2double(get(hObject,'String')) returns contents of SegmentationConfidence as a double
handles.SlideThreshold.Value = str2num(handles.SegmentationConfidence.String);


%% --- Executes when figMain is resized.
function figMain_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figMain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

S.fh = handles.figMain;
% Fix the central figure axis in pixel coordinates
S.ax = gca;
set(S.ax,'unit','pix','position',[(handles.figMain.Position(3)-512)/2 handles.figMain.Position(4)-512 512 512]);


%% --- Executes on button press in ViewerMode.
function ViewerMode_Callback(hObject, eventdata, handles)
% hObject    handle to ViewerMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ViewerMode
if (handles.ViewerMode.Value)
    handles.ViewTracking.Visible = 'On';
    handles.ViewTracking.Enable = 'On';
else
    handles.ViewTracking.Visible = 'Off';
    handles.ViewTracking.Enable = 'Off';
end


%% --- Executes on button press in ViewTracking.
function ViewTracking_Callback(hObject, eventdata, handles)
% hObject    handle to ViewTracking (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hMainGui = getappdata(0, 'hMainGui');
data_dir = getappdata(hMainGui, 'data_dir');
if(isempty(getappdata(hMainGui, 'trajectory')))
    setappdata(hMainGui, 'CoM', []);
    setappdata(hMainGui, 'trajectory', []);
    setappdata(hMainGui, 'norm_trajectory', []);
    setappdata(hMainGui, 'start_frame', []);
end
CoM = getappdata(hMainGui, 'CoM');
trajectory = getappdata(hMainGui, 'trajectory');
norm_trajectory = getappdata(hMainGui, 'norm_trajectory');

pos_bs = strfind(data_dir,'Data');
sub_dir = data_dir(pos_bs(end)+length('Data'):length(data_dir));
seg_dir = ['./Results/SegmentedImages' sub_dir '/'];
data_dir = ['./Data' sub_dir '/'];
if(~isempty(dir([data_dir '*.tif'])))
    img_list = dir([data_dir '*.tif']);
else
    img_list = dir([data_dir '*.bmp']);
end
output_dir = ['./Results/Tracking/' sub_dir '/'];
if(~exist(output_dir))
    mkdir(output_dir);
end

if(length(handles.ChooseDisplay.String) <3)
    handles.ChooseDisplay.String{3} = 'Tracking';
end

nLegs = 8 - handles.L4.Value - handles.R4.Value;
if (nLegs > 6)
    colors = [1 0 0; 0 1 0; 0 0 1; 1 0 0.5; 1 1 0; 0 1 1; 1 0 1; 1 0.5 0];
    legs_id = {'L1', 'L2', 'L3','L4', 'R1', 'R2', 'R3','R4'};
else
    colors = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 0 1 1; 1 0 1; 1 0 0.5; 1 0.5 0];
    legs_id = {'L1', 'L2', 'L3', 'R1', 'R2', 'R3'};
end

end_frame = length(dir([seg_dir 'img_*.png']));

if (strcmp(handles.ViewTracking.String, 'View Tracking'))
    handles.ViewTracking.String = 'Pause Viewing';
else
    handles.ViewTracking.String = 'View Tracking';
end
while (strcmp(handles.ViewTracking.String, 'Pause Viewing'))
    setappdata(hMainGui, 'current_frame', getappdata(hMainGui, 'current_frame') + 1);
    i = getappdata(hMainGui, 'current_frame');
    if (i>end_frame)
        setappdata(hMainGui, 'current_frame', end_frame);
        break;
    end
    if(~isempty(trajectory))
        handles.text.String = ['Current frame: ' num2str(i)];
        handles.GoToDropdown.Value = i;
        
        handles.Tracking_Text.String = '';
        imshow([data_dir img_list(i).name]);
        hold on;
        for kk = 1:nLegs
            if(eval(['handles.' legs_id{kk} '.Value']))
                handles.Tracking_Text.String = [handles.Tracking_Text.String; sprintf('Leg %2s:        ', legs_id{kk})];
            end
            if(trajectory(i,kk,1)~=0)
                handles.Tracking_Text.String = [handles.Tracking_Text.String; sprintf('Leg %2s: %3d %3d', legs_id{kk}, trajectory(i,kk,1), trajectory(i,kk,2))]; 
                scatter(trajectory(i,kk,1),trajectory(i,kk,2),'w');
                text(trajectory(i,kk,1),trajectory(i,kk,2),legs_id{kk},'Color',[colors(kk,1) colors(kk,2) colors(kk,3)],'FontSize',14);
            else
                handles.Tracking_Text.String = [handles.Tracking_Text.String; sprintf('Leg %2s:        ', legs_id{kk})];
            end
        end
        hold off;
        drawnow
        
        % Update Display
        nLegs = 8 - handles.L1.Value - handles.L2.Value - handles.L3.Value - handles.R1.Value - handles.R2.Value - handles.R3.Value - handles.L4.Value - handles.R4.Value;
        if(length(str2num(handles.Tracking_Text.String(:,9:11))) < nLegs)
            set(handles.Missing, 'BackgroundColor', [1 0 0]); 
        else
            set(handles.Missing, 'BackgroundColor', [0 1 0]); 
        end
    end
end


%% --- Executes on button press in DataProcess.
function DataProcess_Callback(hObject, eventdata, handles)
% hObject    handle to DataProcess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hMainGui = getappdata(0, 'hMainGui');

if(isempty(getappdata(hMainGui, 'data_dir')))
    set(handles.Text_Display, 'String', 'Analysing tracked data.');
    pause(0.5);
    tracking_dir = uigetdir('./Results/Tracking/');
    
    defaultfps = {'1000'};
    fps = inputdlg({'Recording FPS for dataset?'},'Input Recording FPS',1,defaultfps);
    fps = str2num(fps{:});

    defaultscale = {'11.0'};
    scale = inputdlg({'Image Size (mm)?'},'Image Size Input',1,defaultscale);
    scale = str2double(scale{:});

    DataProcessing_New(fps,tracking_dir,scale);
    set(handles.Text_Display, 'String', 'Data Processing complete.');
else
    data_dir = getappdata(hMainGui, 'data_dir');
    
    defaultfps = {'1000'};
    fps = inputdlg({'Recording FPS for dataset?'},'Input Recording FPS',1,defaultfps);
    fps = str2num(fps{:});

    defaultscale = {'11.0'};
    scale = inputdlg({'Image Size (mm)?'},'Image Size Input',1,defaultscale);
    scale = str2double(scale{:});

    pos_bs = strfind(data_dir,'Data');
    sub_dir = data_dir(pos_bs(end)+length('Data'):length(data_dir));
    tracking_dir = ['./Results/Tracking' sub_dir '/'];
    if(~exist([tracking_dir 'CoM.mat']))
        handles.Text_Display.String = 'Please ensure that tracking has been completed first.';
    else
        set(handles.Text_Display, 'String', 'Analysing tracked data.');
        pause(0.5);
        DataProcessing_New(fps,tracking_dir,scale);
        set(handles.Text_Display, 'String', 'Data Processing complete.');
    end
end


%% --- Executes on button press in BatchDataProcess.
function BatchDataProcess_Callback(hObject, eventdata, handles)
% hObject    handle to BatchDataProcess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data_dir = uipickfiles('filterspec',[pwd '/Results/Tracking']);
defaultfps = {'1000'};
fps = inputdlg({'Recording FPS for dataset?'},'Input Recording FPS',1,defaultfps);
fps = str2num(fps{:});

defaultscale = {'11.0'};
scale = inputdlg({'Image Size (mm)?'},'Image Size Input',1,defaultscale);
scale = str2double(scale{:});

try
    if(data_dir==0)
        return;
    end
catch
    
end
if(isempty(data_dir))
    return;
end
for i = 1 : length(data_dir)
    AnalyseFolder(data_dir{i},fps,scale);
end

function AnalyseFolder(data_dir,fps,scale)
    pos_bs = strfind(data_dir,'Results');
    sub_dir = data_dir(pos_bs(end)+length('Results/Tracking/'):length(data_dir));
    output_dir = ['./Results/Tracking/' sub_dir '/'];
    if(~exist([output_dir 'trajectory.mat']))
        dir_fn = dir(output_dir);
        for i_dir = 1:length(dir_fn)-2
            AnalyseFolder([output_dir '/' dir_fn(i_dir+2).name],fps,scale);
            pause(1);
        end
    else
        try    
            DataProcessing_New(fps,output_dir,scale);
        catch MExc
            fprintf('An error occured during data processing for the folder: %s.\n', data_dir);
            disp('Execution will continue.');
            disp(MExc.message);
        end
    end


%% --- Executes on key press with focus on figMain or any of its controls.
function figMain_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figMain (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
if(strcmp(handles.nextframe.Enable,'On'))
    switch eventdata.Key
        case 'leftarrow', prevframe_Callback(hObject, eventdata, handles);
        case 'rightarrow', nextframe_Callback(hObject, eventdata, handles);
        case 'uparrow', prevframe_Callback(hObject, eventdata, handles);
        case 'downarrow', nextframe_Callback(hObject, eventdata, handles);
    end
end


%% --- Executes when user attempts to close figMain.
function figMain_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figMain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Construct a questdlg with three options
hMainGui = getappdata(0, 'hMainGui');

choice = questdlg('Are you sure?', ...
    'FLLIT', ...
    'Save and Exit', 'Exit without Saving', 'Cancel', 'Cancel');

% Handle response
switch choice
    case 'Save and Exit'
        if(~isempty(getappdata(hMainGui, 'trajectory')))        
            data_dir = getappdata(hMainGui, 'data_dir');
            pos_bs = strfind(data_dir,'Data');
            sub_dir = data_dir(pos_bs(end)+length('Data'):length(data_dir));
            output_dir = ['./Results/SegmentedImages' sub_dir '/'];
            trajectory = getappdata(hMainGui, 'trajectory');
            norm_trajectory = getappdata(hMainGui, 'norm_trajectory');
            CoM = getappdata(hMainGui, 'CoM');
            save([output_dir 'CoM.mat'],'CoM');
            save([output_dir 'trajectory.mat'],'trajectory');
            save([output_dir 'norm_trajectory.mat'],'norm_trajectory');
            csvwrite([output_dir 'CoM.csv'],CoM);
            csvwrite([output_dir 'trajectory.csv'],trajectory);
            csvwrite([output_dir 'norm_trajectory.csv'],norm_trajectory);
            display('Saving... Exiting...');
        else
            display('Nothing to save! Exiting...');
        end
        T = timerfind;
        if (~isempty(T))
            stop(T);
            delete(T);
        end
        clear T
        delete(hObject);
    case 'Exit without Saving'
        display('Exiting...');
        T = timerfind;
        if (~isempty(T))
            stop(T);
            delete(T);
        end
        clear T
        delete(hObject);
    case 'Cancel'
        return;
end


% --- Executes on slider movement.
function ImageSizeSlider_Callback(hObject, eventdata, handles)
% hObject    handle to ImageSizeSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.ImageSize.String = num2str(hObject.Value);


function ImageSize_Callback(hObject, eventdata, handles)
% hObject    handle to ImageSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ImageSize as text
%        str2double(get(hObject,'String')) returns contents of ImageSize as a double
handles.ImageSizeSlider.Value = str2num(hObject.String);


% --- Executes on button press in ResetImageSize.
function ResetImageSize_Callback(hObject, eventdata, handles)
% hObject    handle to ResetImageSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
bodylength = str2num(handles.BodyLength.String);
if(bodylength>0)
    handles.ImageSizeSlider.Value = round(128/bodylength*11*2)/2;
    handles.ImageSize.String = num2str(round(128/bodylength*11*2)/2);
else
    handles.ImageSizeSlider.Value = 11;
    handles.ImageSize.String = num2str(11);
end
