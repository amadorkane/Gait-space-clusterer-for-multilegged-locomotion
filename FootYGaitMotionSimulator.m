%% compute a set of gait phase difference data from one or more model gaits
% add random noise to each model gait to create a perturbed set of data

% initialization
clear all; close all; clc;

%% get phase difference spreadsheet
[gaitdatafile, filepath] = uigetfile('*.*','Select file with model gaits');
cd(filepath);
% get the corefilename to use in naming other files
[~,corefilename,~] = fileparts(gaitdatafile);
% create full filename
datafilename = fullfile(filepath,gaitdatafile);

%% get filename and data from the master spreadsheet file
inputdata = readcell(datafilename);

%% select model to simulate
header = inputdata(1,:);
model_list = inputdata(2:end,1);  % model names
nmodels = numel(model_list);          % number models
% model phase differences
model_data = cell2mat(inputdata(2:end,2:end));
Nlegs = size(model_data,2) + 1;

[model_idx,~] = listdlg("PromptString",'select gait to simulate',...
    "SelectionMode","single",...
    "InitialValue",1,"ListString",model_list);

% get phase differences and name of the selected model
model_data = model_data(model_idx,:);
model_name = model_list{model_idx};

outputfile = fullfile(filepath,[model_name,num2str(Nlegs),'leg_data_xypts.csv']);

%% select simulation parameters

% get number of legs
answer = inputdlg({'stride period(frames)','number strides'},...
    'simulation parameters',[1 40],{'50','5','2.9'});

% number of random samples of each perturbed model gait
stride_period = str2num(answer{1}); 
nstrides = str2num(answer{2}); 

nframes = nstrides*stride_period;  % total number of frames

%% create perturbed simulated idealized model gaits

simdata = [];

% generate array of all phase differences
for j = 1:nframes
    for i = 1:Nlegs - 1
        simphasediff(j,i) = model_data(i);
    end
end

% compute cumulative phase differences along rows to get the phase
% difference of each of the i+1th foot relative to i=1 foot
cumsimphasediff = cumsum(simphasediff,2);

%% compute y motion for each foot
FootXY = zeros(nframes,2*Nlegs);
for j = 1:nframes
    for i = 1:Nlegs
        if i == 1
            FootXY(j,2*i) = cos(j*2*pi/stride_period);  % phase = 0 by definition
        else
            FootXY(j,2*i) = cos(j*2*pi/stride_period + 2*pi*cumsimphasediff(j,i-1)); % phase shifted
        end
        FootXY(j,2*i-1) = 0;             % all x = 0
    end
end

%% save the data in the specified format of the Nlegs-1 phase differences
% in the first Nlegs-1 columns and the model gait names in the last column

% add the body COM, Cranial, Caudal end XY to this set for using to test the
% PhaseDiffCalculator.m code

FootXY = [zeros(nframes,6),FootXY];
% create header
header = {};
for i = 1:(Nlegs + 3)
    header = [header ['pt',num2str(i),'_cam1_X'],['pt',num2str(i),'_cam1_Y']];
end

outputarray = [header; num2cell(FootXY)];

writecell(outputarray,outputfile);

%% plot the foot motion along y 
% (ignoring first 6 columns = body COM, Cranial, Caudal end XY)
for i = 1:Nlegs
    legnames{i} = ['leg',num2str(i)];
end
stackedplot(FootXY(:,8:2:end),"DisplayLabels",legnames);
xlabel('frames'); title(model_name);

%% end of main program

%% ************************************************************************