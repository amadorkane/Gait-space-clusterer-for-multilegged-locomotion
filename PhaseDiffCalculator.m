%% PhaseDiffCalculator.m
% Reads in tracked coordinates of the body COM, head, caudal end, and 
% tracked distal leg endpoints in DLTdv format and uses them to compute the phase
% differences between adjacent legs
%
% Input file: 
% (corefilename)_data_xypts.csv data file in DLTdv format
%
% MUST HAVE the body cranial-caudal axis to be oriented along the image
% vertical axis, where the caudal-to-cranial direction points upward
%
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
%
% Output: 
% (corefilename)_phasediff.csv file with phase differences between Nlegs-1
% adjacent leg pairs for clustering and other gait analysis

%% Initialization
close all; clear all; clc;
polyorder = 2;              % use quadratic polynomial in computing slopes
%---------------------
% video & data analysis
framerate = 500;            % video framerate = sampling frequency

tsmooth = 25;                    % smoothing by 25 ms (40 Hz, 12.5 frames)
% using a quadratic running fit
% reduce tsmooth if stride frequencies > 20 Hz
detrend_polyorder = 2;      % polynomial order for detrending data: 0 = constant
% 1 = linear; 2 = quadratic, etc.

%% read in track data

[trackfile,path] = uigetfile('*_data_xypts.csv','select foot coordinate file');

trackfilename = fullfile(path,trackfile);

% create short file name to use in creating output files
[~,corefilename,~] = fileparts(trackfilename);
corefilename((end-10):end) = [];

% create output file for phase differences
PhaseDiffFilename = fullfile(path,[corefilename,'_phasediff.csv']);
if isfile(PhaseDiffFilename)             % delete if it already exists
    delete(PhaseDiffFilename);
end

%% read in the track coordinates
track_coords = readcell(trackfilename);

% parse into specific variables
% remove header (row 1)
track_coords(1,:) = [];
% convert to numerical
track_coords = cell2mat(track_coords);
%% extract tracks for each tracked landmark
% for computing phase, COM needs to be fixed at the image center in the 
% body-fixed frame;  use to correct foot x,y if necessary later
COM = track_coords(:,1:2);
% Head = track_coords(:,3:4);    % not needed
% Caudal = track_coords(:,5:6);  % not needed
nframes = size(track_coords,1);
% find all rows that are all NaN so we can ignore:
idx_nan = all(isnan(track_coords),2);
diff_idx_nan = diff(idx_nan);
numNaNstart = find(diff_idx_nan == -1);
if isempty(numNaNstart) % no all NaN rows at start
    firstframe = 1;
else
    firstframe = numNaNstart(1) + 1; % in case more than one sequence
end
numNaNend = find(diff_idx_nan == 1);
if isempty(numNaNend) % no all NaN rows at start
    lastframe = nframes - firstframe + 1;
else
    lastframe = numNaNend(end); % in case more than one sequence
end
%% input analysis parameters
prompt = {'leg motion smoothing window(ms)',...
    'framerate(frame/s)', 'detrend polynomial order',...
    'tracked +y points up/down (1=Y,0=N)'};
dlgtitle = 'Analysis parameters';
fieldsize = [1 40];
definput = {num2str(tsmooth),num2str(framerate),num2str(detrend_polyorder),'1'};

answer = inputdlg(prompt,dlgtitle,fieldsize,definput);

tsmooth             = str2num(answer{1});
framerate           = str2num(answer{2});
detrend_polyorder   = str2num(answer{3});
switch_y_dir       = ~str2num(answer{4});

dt = 1/framerate;                           % time step
smoothwindow = ceil(tsmooth*1e-3/dt);       % smoothing window for foot motion

%% process foot motions

% get original foot positions in the body fixed frame:
FootXY = track_coords(:,7:end);
% number of legs = number of xy columns/2
Nlegs = round(size(FootXY,2))/2;
% ------------------
% fill in any missing NaN values with a spline
for i = 1:size(FootXY,2)
    FootXY(:,i) = fillmissing(FootXY(:,i),'spline');
end
% ------------------
% save tarsus coords in nicer format:
clear temp;
for j = 1:Nlegs
    temp(:,:,j) = FootXY(:,2*(j-1)+[1,2]);
end
FootXY = temp;                  % dim = (frame,(x,y),tarsus)
% ------------------
% get coordinates relative to body COM
for j = 1:Nlegs
    FootXY(:,:,j) = FootXY(:,2*(j-1)+[1,2]) - COM;
end
% ------------------
if switch_y_dir
    % reverse direction of y coordinates so +y-axis points upward (original
    % image convention has opposite orientation)
    FootXY(:,2,:) = -FootXY(:,2,:);
end
% ------------------
% FootXY now holds the coordinates of each foot in the body-fixed
% coordinate frame where +y = cranial-caudal axis and the origin = COM

% for use in computing oscillation phase below , save the y position of the
% foot along the cranial-caudal axis in the body fixed frame
Footy_cc(:,:) = FootXY(:,2,:);

% Because the legs can be held at somewhat different mean positions for
% different motions, we first detrend the y values; this allows for correct
% normalization of y when computing phases; this removes a long-time value
% of speed that isn't relevant for phase calculations.
Footy_cc(firstframe:lastframe,j) = detrend(Footy_cc(firstframe:lastframe,j),detrend_polyorder);

%% compute velocity along the +y cranial-caudal direction for each leg
% using quadratic polynomial local fits
for j = 1:Nlegs
    FootvY_cc(firstframe:lastframe,j) = ...
        movingslope(FootXY(firstframe:lastframe,2,j),smoothwindow,polyorder,dt);
end

%% plot y and vy vs time
keepgoing = 1;
while keepgoing
    figure('Position',[10 50 1400 900]);
    tiledlayout(Nlegs,2,'TileSpacing','compact','Padding','compact');
    for j = 1:Nlegs
        nexttile;
        plot(firstframe:lastframe,normalize(Footy_cc(firstframe:lastframe,j)));
        hold on;
        ylabel(['y leg ',num2str(j)]); 
        plot([firstframe,lastframe],[0,0],'k-');
        nexttile;
        plot(firstframe:lastframe,normalize(FootvY_cc(firstframe:lastframe,j)));
        hold on;
        ylabel(['v_y leg ',num2str(j)]); 
        plot([firstframe,lastframe],[0,0],'k-');
    end
    xlabel('frame');
    % select whether or not to compute using a limited range of the input
    % data -- phase differences only make sense if you have an oscillatory
    % signal, because in the next step we're going to use a normalized
    % signal
    answer = questdlg('The y and v_y data must be oscillatory for phase difference analysis--select new first and last frames?', ...
    	'Choices', ...
    	'Yes','No','No');
    switch answer
        case 'Yes'
            prompt = {'first frame:','last frame:'};
            dlgtitle = 'Data must be oscillatory for phase differences to make sense';
            fieldsize = [1 40];
            definput = {num2str(firstframe),num2str(lastframe)};

            answer = inputdlg(prompt,dlgtitle,fieldsize,definput);
            % only compute phase differences for these frames
            firstframe      = str2num(answer{1});
            lastframe       = str2num(answer{2});
            clf;
        otherwise
            keepgoing = 0;
    end
end
%% compute, analyze & plot oscillation phase for each leg

for j = 1:Nlegs
    % Compute phase angle from standardized values of the position and
    % velocity in the comoving frame where each canonical variable is
    % computed along the direction of COM velocity at each instant of time;
    % Ref: see Revzen et al., 2013 equation 1 & surrounding text. 
    % 
    % Compute the oscillation phase from equations:
    % y = A cos(wt); v_y = -Awsin(wt)
    % So we need to standardize both y and v_y (analogous to dividing by A 
    % and Aw) and include a factor of -1: 
    % tan wt = sin wt/cos wt = v_y/(-Aw)/y/A 
    % so the phase is:
    % wt = atan(v_x/(-Aw)/x/A) 
    % 
    % Here we use normalize so y and vy have zero mean, unit stdev 

    osc_phase(:,j) = atan2(-normalize(FootvY_cc(firstframe:lastframe,j)),...
        normalize(Footy_cc(firstframe:lastframe,j)));

    % the range for atan2 0 to +pi, then -pi to 0, but we need all positive
    % values [0, 2pi]
    %
    % if osc_phase (rad) < 0, then add 2pi so range = [0,2pi]
    % map onto 0 to 2pi by adding 2
    % pi to negative values
    osc_phase((osc_phase(:,j) < 0),j) = 2*pi + osc_phase((osc_phase(:,j)<0),j);

    % unwrap phases (rad) so each phase varies continuously (resolves jumps
    % with magnitude > pi by adding or subtracting 2*pi until the jump is
    % smaller than pi in magnitude

    osc_phase(:,j) = unwrap(osc_phase(:,j));

    % divide by 2pi to convert from radians to cycles 

    osc_phase(:,j) = osc_phase(:,j)/(2*pi); % cycles
end

% compute the phase differences osc_phase_diff in units of cycles:

for j = 1:(Nlegs-1)
    osc_phase_diff(:,j) = osc_phase(:,j+1) - osc_phase(:,j);
end

% Restrict the phase difference range to [0,1] cycle = [0,2*pi]
%
% if phi < 0, a phase difference of phi (in cycles) is the same as 1 + phi
% cycles; if the phase difference is > 1, then subtract 1 from it
for j = 1:(Nlegs-1)
    temp = osc_phase_diff(:,j);
    osc_phase_diff(:,j) = mod(temp,1);
end

%% save oscillation phase differences between adjoining intact legs 

% create header
header = {};
for j = 1:(Nlegs-1)
    header = [header ['phase diff ',num2str(j),'(cycle)']];
end
header = [header 'frame'];

outputarray = num2cell([osc_phase_diff (firstframe:lastframe)']);
outputarray = [header; outputarray];

% write phase differences to file
writecell(outputarray,PhaseDiffFilename);

%% END of main program

%% ************************************************************************
% supporting functions
%% ************************************************************************

function Dvec = movingslope(vec,supportlength,modelorder,dt)
% movingslope: estimate local slope for a sequence of points, using a sliding window
% usage: Dvec = movingslope(vec)
% usage: Dvec = movingslope(vec,supportlength)
% usage: Dvec = movingslope(vec,supportlength,modelorder)
% usage: Dvec = movingslope(vec,supportlength,modelorder,dt)
%
%
% movingslope uses filter to determine the slope of a curve stored
% as an equally (unit) spaced sequence of points. A patch is applied
% at each end where filter will have problems. A non-unit spacing
% can be supplied.
%
% Note that with a 3 point window and equally spaced data sequence,
% this code should be similar to gradient. However, with wider
% windows this tool will be more robust to noisy data sequences.
%
% arguments: (input)
%  vec - row of column vector, to be differentiated. vec must be of
%        length at least 2
%
%  supportlength - (OPTIONAL) scalar integer - defines the number of
%        points used for the moving window. supportlength may be no
%        more than the length of vec.
%
%        supportlength must be at least 2, but no more than length(vec)
%
%        If supportlength is an odd number, then the sliding window
%        will be central. If it is an even number, then the window
%        will be slid backwards by one element. Thus a 2 point window
%        will result in a backwards differences used, except at the
%        very first point, where a forward difference will be used.
%
%        DEFAULT: supportlength = 3
%
%  modelorder - (OPTIONAL) - scalar - Defines the order of the windowed
%        model used to estimate the slope. When model order is 1, the
%        model is a linear one. If modelorder is less than supportlength-1.
%        then the sliding window will be a regression one. If modelorder
%        is equal to supportlength-1, then the window will result in a
%        sliding Lagrange interpolant.
%
%        modelorder must be at least 1, but not exceeding
%        min(10,supportlength-1)
%
%        DEFAULT: modelorder = 1
%
%  dt - (OPTIONAL) - scalar - spacing for sequences which do not have
%        a unit spacing.
%
%        DEFAULT: dt = 1
%
% arguments: (output)
%  Dvec = vector of derivative estimates, Dvec will be of the same size
%        and shape as is vec.
% 
%
% Example:
%  Estimate the first derivative using a 7 point window with first through
%  fourth order models in the sliding window. Note that the higher order
%  approximations provide better accuracy on this curve with no noise.
%  
%  t = 0:.1:1;
%  vec = exp(t);
%
%  Dvec = movingslope(vec,7,1,.1)
%  Dvec =
%  Columns 1 through 7
%    1.3657  1.3657  1.3657  1.3657  1.5093  1.668  1.8435
%  Columns 8 through 11
%    2.0373  2.0373  2.0373  2.0373
%
%  Dvec = movingslope(vec,7,2,.1)
%  Dvec =
%  Columns 1 through 7
%    0.95747 1.0935  1.2296  1.3657  1.5093  1.668  1.8435
%  Columns 8 through 11
%    2.0373  2.2403  2.4433  2.6463
%
%  Dvec = movingslope(vec,7,3,.1)
%  Dvec =
%  Columns 1 through 7
%    1.0027  1.1049  1.2206  1.3498  1.4918  1.6487  1.8221
%  Columns 8 through 11
%    2.0137  2.2268  2.4602  2.7138
%
%  Dvec = movingslope(vec,7,4,.1)
%  Dvec =
%    Columns 1 through 7
%    0.99988 1.1052  1.2214  1.3498  1.4918  1.6487  1.8221
%  Columns 8 through 11
%    2.0137  2.2255  2.4597  2.7181
%
%
% Example:
%  Estimate the slope of a noisy curve, using a locally quadratic
%  approximation. In this case, use a straight line so that we know
%  the true slope should be 1. Use a wide window, since we have
%  noisy data.
%  
%  t = 0:100;
%  vec = t + randn(size(t));
%  Dvec = movingslope(vec,10,2,1)
%  mean(Dvec)
%  ans = 
%     1.0013
%  std(Dvec)
%  ans =
%     0.10598
%
%  By way of comparison, gradient gives a much noisier estimate
%  of the slope of this curve.
%
%  std(gradient(vec))
%  ans =
%     0.69847
%
%
% Example:
%  As a time test, generate random data vector of length 500000.
%  Compute the slopes using a window of width 10.
%
%  vec = rand(1,500000);
%  tic
%  Dvec = movingslope(vec,10,2);
%  toc
%
%  Elapsed time is 0.626021 seconds.
%
%
% See also: gradient
%
%
% Author: John D'Errico
% e-mail: woodchips@rochester.rr.com
% Release: 1.0
% Release date: 10/19/07
% Cite as:  John D'Errico (2021). Movingslope
%  (https://www.mathworks.com/matlabcentral/fileexchange/16997-movingslope),
%  MATLAB Central File Exchange. Retrieved August 11, 2021.
%**************************************************************************
% how long is vec? is it a vector?
if (nargin==0)
  help movingslope
  return
end
if ~isvector(vec)
  error('vec must be a row or column vector')
end
n = length(vec);                    % n = number of input data points

% -------------------------------------------------------------------------
% supply defaults
if (nargin<4) || isempty(dt)
  dt = 1;
end
if (nargin<3) || isempty(modelorder)
  modelorder = 1;
end
if (nargin<2) || isempty(supportlength)
  supportlength = 3;
end
% check the parameters for problems
if (length(supportlength)~=1) || (supportlength<=1) || (supportlength>n) || (supportlength~=floor(supportlength))
  error('supportlength must be a scalar integer, >= 2, and no more than length(vec)')
end
if (length(modelorder)~=1) || (modelorder<1) || (modelorder>min(10,supportlength-1)) || (modelorder~=floor(modelorder))
  error('modelorder must be a scalar integer, >= 1, and no more than min(10,supportlength-1)')
end
if (length(dt)~=1) || (dt<0)
  error('dt must be a positive scalar numeric variable')
end

% -------------------------------------------------------------------------
% now build the filter coefficients to estimate the slope
if mod(supportlength,2) == 1
  parity = 1; % odd parity
else
  parity = 0;
end
s = (supportlength-parity)/2;       % s = approx. half window size 
t = ((-s+1-parity):s)';             % t = window for filter calculations
coef = getcoef(t,supportlength,modelorder);

% -------------------------------------------------------------------------
% Apply the filter to the entire vector
f = filter(-coef,1,vec);
Dvec = zeros(size(vec));
Dvec(s+(1:(n-supportlength+1))) = f(supportlength:end);
% --------------------------
% patch each end
vec = vec(:);
for i = 1:s
  % patch the first few points (1 to s) so they only involve fits to the
  % the appropriate points
  t = (1:supportlength)' - i;
  coef = getcoef(t,supportlength,modelorder);
  
  Dvec(i) = coef*vec(1:supportlength);
  
  % patch the end points ((n - s + 1) to n) so they only involve fits over
  % the appropriate points
  if i<(s + parity)
    t = (1:supportlength)' - supportlength + i - 1;
    coef = getcoef(t,supportlength,modelorder);
    Dvec(n - i + 1) = coef*vec(n + (0:(supportlength-1)) + 1 - supportlength);
  end
end

% -------------------------------------------------------------------------
% scale by the supplied spacing
Dvec = Dvec/dt;
% all done
end % mainline end
% =========================================================
% subfunction, used to compute the filter coefficients
function coef = getcoef(t,supportlength,modelorder)
    % Note: bsxfun would have worked here as well, but some people
    % might not yet have that release of matlab.
    A = repmat(t,1,modelorder+1).^repmat(0:modelorder,supportlength,1);
    pinvA = pinv(A);
    % we only need the linear term
    coef = pinvA(2,:);
end % nested function end

%% ************************************************************************