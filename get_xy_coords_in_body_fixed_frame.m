%% get_xy_coords_in_body_fixed_frame.m
% converts DLTdv tracked x,y coordinates into foot (i.e., leg distal end)
% positions in a body center of mass (COM) fixed frame with +y (up on the
% image) along the body cranial-caudal axis
% 
% Input spreadsheet file must be formatted as a DLTdv
% (corefilename)_data_xypts.csv spreadsheet
% row 1 = header (ignored)
% order of tracked points: 
% pt1 = body COM (col 1-2)
% pt2 = head (cranial end) (col 3-4) (not used)
% pt3 = caudal end (col 5-6)  (not used)
% pt4, ... = tracked distal leg endpoints (col 7 - end), starting with the 
% foremost left leg and going counterclockwise (in dorsal view) to the 
% foremost right leg
%
% Input parameters:
% > smoothing window in millisec for computing time derivaties as slopes
% flag to choose whether +y is up or down on the image, depending on
% tracking convention (program can flip to correct)
% > framerate of original video
% > whether the upward direction on the original video) corresponds to +y
% or -y (different image analysis conventions)

% framerate
%
% Output: spreadsheet for inputting to PhaseDiffCalculator.m to compute
% phase differences between adjacent legs
% (corefilename)_xypts.csv file with phase differences between Nlegs-1
% adjacent leg pairs for clustering and other gait analysis

%% Initialization
clc; clear all; close all

markertype = {'o' 's' 'd' '^' 'v' '*' 'h' 'p' '+'};

%---------------------
% video & data analysis
framerate = 500;            % video framerate = sampling frequency
tsmooth = 50;               % smooth body coords by 1/20Hz=50 ms=25 frames

%% read in track data

[trackfile,path] = uigetfile('*_data_xypts.csv','select foot coordinate file');

trackfilename = fullfile(path,trackfile);

% create short file name to use in creating output files
[~,corefilename,~] = fileparts(trackfilename);
corefilename((end-10):end) = [];

% create output xy coordinate file
outputfilename = fullfile(path,[corefilename,'_fix_data_xypts.csv']);
if isfile(outputfilename)             % delete if it already exists
    delete(outputfilename);
end

%% input analysis parameters
prompt = {'leg motion smoothing window(ms)',...
    'framerate(frame/s)', ...
    'tracked +y points up/down (1=Y,0=N)'};
dlgtitle = 'Analysis parameters';
fieldsize = [1 40];
definput = {num2str(tsmooth),num2str(framerate),'1','0'};

answer = inputdlg(prompt,dlgtitle,fieldsize,definput);

tsmooth         = str2num(answer{1});
framerate       = str2num(answer{2});
switch_y_dir    = ~str2num(answer{3});

dt = 1/framerate;                           % time step
smoothwindow = ceil(tsmooth*1e-3/dt);       % smoothing window for foot motion

%% read in the track coordinates
track_coords = readcell(trackfilename);

% parse into specific variables
% save, then remove header (row 1)
header = track_coords(1,:);
track_coords(1,:) = [];
track_coords = cell2mat(track_coords);  % convert from cell array to matrix

% ------------------
if switch_y_dir
    % reverse direction of y coordinates (even columns) so +y-axis points 
    % upward (some image conventions have opposite orientation)
    track_coords(:,2:2:end) = -track_coords(:,2:2:end);
end
%
ncols = size(track_coords,2);        % number tracked x,y values
nframes = size(track_coords,1);      % number frames

%% find & remove any rows with any of these conditions:
% >= half NaN's (i.e. not tracked) at the start 
% and end of the dataset
idx = find(sum(isnan(track_coords),2) >= ncols/3);   % get indices of rows with >= half NaN's 
% trim off untracked rows at start
trip_start_untracked_rows = 1;
while trip_start_untracked_rows
    if ismember(1,idx)                  % first frame is all NaN
        track_coords(1,:) = [];         % delete first row
        idx = find(sum(isnan(track_coords),2) >= ncols/3);   % get indices of rows with >= half NaN's
    else
        trip_start_untracked_rows = 0;  % done trimming
    end
end
nframes = size(track_coords,1);         % number tracked frames
% trim off untracked rows at end
trip_end_untracked_rows = 1;
while trip_end_untracked_rows
    if ismember(nframes,idx)              % last frame is all NaN
        track_coords(end,:) = [];         % delete first row
        idx = find(sum(isnan(track_coords),2) >= ncols/3);   % get indices of rows with >= half NaN's
        nframes = size(track_coords,1);
    else
        trip_end_untracked_rows = 0;  % done trimming
    end
end
nframes = size(track_coords,1);         % number tracked frames

%% extract tracks for each tracked landmark
% body center-of-mass
COM = track_coords(:,1:2);
Head = track_coords(:,3:4);    
Caudal = track_coords(:,5:6);

% get foot (leg distal end) coordinates
FootXY = track_coords(:,7:end);
% number of legs = number of xy columns/2
Nlegs = round(size(FootXY,2))/2;

% ------------------
% save foot coords in nicer format
% and subtract off COM so it's at the origin
clear temp;
for j = 1:Nlegs
    temp(:,:,j) = FootXY(:,2*(j-1)+[1,2]);
end
FootXY = temp;                  % dim = (frame,(x,y),leg)

% fill in missing values
COM     = fillmissing(COM,'movmedian',smoothwindow);
Head    = fillmissing(Head,'movmedian',smoothwindow); 
Caudal  = fillmissing(Caudal,'movmedian',smoothwindow);

% smooth data
COM = smoothdata(COM,smoothwindow,'rloess',"omitmissing");
Head = smoothdata(Head,smoothwindow,'rloess',"omitmissing");
Caudal = smoothdata(Caudal,smoothwindow,'rloess',"omitmissing");

%% compute normalized cranial caudal (cc) axis
ccaxis = (Head - Caudal); ccaxis = ccaxis./vecnorm(ccaxis);

% compute orientation of the cranial-caudal axis (yaw)
for i = 1:nframes
    yaw_angle(i) = atan2d(ccaxis(i,2),ccaxis(i,1));  % yaw in deg CCW from +x
end

% unwrap before smoothing or else it might smooth over 0 to > or < 360 deg
% jumps in the data
yaw_angle = unwrap(yaw_angle*(pi/180))*(180/pi);
yaw_angle = smooth(yaw_angle,smoothwindow,'rloess',"omitmissing");

for i = 1:nframes
    % we want to rotate all coordinates so the cc axis agrees with +y
    % this requires rotating all coordinates by the following angles
    % rotm = rotation matrix to rotate xy coords into the correct orientation
    if yaw_angle(i) <= 0      % correct negative yaw
        yaw_angle(i) = yaw_angle(i)+360;
    end
    % now yaw angle is CCW from the +x axis
    R = rotz(90-yaw_angle(i));
    % remove the 3rd (z) row and column so it's a 2D rotation matrix
    R(:,3) = []; R(3,:) = [];
    rotm{i} = R;
end
display('done computing rotational matrix');
%% process foot motions

% subtract off COM so it's at the origin
FootXY = FootXY - COM;
Head = Head - COM;  Caudal  = Caudal - COM;
COM = COM - COM;

display('have coordinates in body-fixed frame')
%% get x,y coordinates in frame where the cranial-caudal axis points along +y
% and the COM is the origin

% Rotate each coordinate
display('rotating coordinates');
for i = 1:nframes
    temp = rotm{i}*Head(i,:)';
    Headrot(i,:) = temp';
    temp =  rotm{i}*Caudal(i,:)';
    Caudalrot(i,:) = temp';
    for j = 1:Nlegs
        XY = squeeze(FootXY(i,:,j))';
        XYrot = rotm{i}*XY;
        FootXYrot(i,:,j) = XYrot';
    end
end

% FootXY now holds the coordinates of each foot in the body-fixed
% coordinate frame where +y = cranial-caudal axis and the origin = COM

%% plot transformed coordinates

figure;
colororder("glow12");  % change plot colors

plot(Head(:,1),Head(:,2),'^'); hold on;
plot(Caudal(:,1),Caudal(:,2),'v');

for j = 1:Nlegs
    plot(FootXY(:,1,j),FootXY(:,2,j),'o');
end
title('original coordinates')
axis equal;

figure;
colororder("glow12");  % change plot colors

plot(Headrot(:,1),Headrot(:,2),'^'); hold on;
plot(Caudalrot(:,1),Caudalrot(:,2),'v');

for j = 1:Nlegs
    plot(FootXYrot(:,1,j),FootXYrot(:,2,j),'o');
end
title('rotated coordinates')
axis equal;

figure;
yaw_angle = unwrap(yaw_angle*(pi/180))*(180/pi);
plot(yaw_angle);
ylabel('yaw (deg)');xlabel('frame')

%% save the coordinates in the new frame
outputarray =  {};

% create single array of all foot variables again
for j = 1:Nlegs
    allFootXYrot(:,2*(j-1)+[1,2]) = FootXYrot(:,:,j);
end
allcoords = [COM Headrot Caudalrot allFootXYrot];
outputarray = [header;num2cell(allcoords)];

writecell(outputarray,outputfilename);

%% end of main program

%% ************************************************************************
% plot y(frame,leg) for each leg vs time
function plot_y_vs_t(y)
    figure('Position',[10 50 1400 900]);
    Nlegs = size(y,2);
    tiledlayout(Nlegs,1,'TileSpacing','compact','Padding','compact');
    for j = 1:Nlegs
        nexttile;
        plot(y(:,j));
        ylabel(['y leg ',num2str(j)]);
    end
    xlabel('frame');
end
%% ************************************************************************
