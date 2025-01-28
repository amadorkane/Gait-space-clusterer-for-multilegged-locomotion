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

outputfile = fullfile(filepath,[corefilename,'_sim.csv']);

%% get filename and data from the master spreadsheet file
inputdata = readcell(datafilename);

%% select simulation parameters

% get number of legs
answer = inputdlg({'number samples each model','von Mise kappa = '},'simulation parameters',[1 40],{'1000','2.9'});

% number of random samples of each perturbed model gait
nsamples = str2num(answer{1});    
% kappa = von Mises distribution width parameter
% measured values for spiders: kappa = 2.9/2.6, 1.6, 4.6 (mean/median, 10%, 90% percentile)
kappa = str2num(answer{2});    

%% input model
header = inputdata(1,:);
model_list = inputdata(2:end,1);  % model names

nmodels = numel(model_list);          % number models

% model phase differences
model_data = cell2mat(inputdata(2:end,2:end));
Nlegs = size(model_data,2) + 1;

%% create perturbed simulated idealized model gaits

simdata = [];

% generate perturbed data
simmodel = {};
for j = 1:nmodels
    for i = 1:Nlegs - 1
        perturbed_model_data(j,:,i) = perturb_gait(model_data(j,i),kappa,nsamples);
    end
    simdata = [simdata;squeeze(perturbed_model_data(j,:,:))];
    modelname    = cell(nsamples,1); modelname(:) = {model_list{j}};
    simmodel = [simmodel;modelname]; %*% make up a list with nsamples copies of hte names
end

%% save the data in the specified format of the Nlegs-1 phase differences
% in the first Nlegs-1 columns and the model gait names in the last column

header_output = [header(1,2:end),header(1,1)];
outputarray = [header_output; num2cell(simdata) simmodel];

writecell(outputarray,outputfile);

fprintf('*********************************************************\n');
%% end of main program

%% ************************************************************************
function alpha = perturb_gait(mu,kappa,nsamples)
    % inputs:  mu = preferred direction of angle in cycles
    % kappa = width parameter in von Mises distribution
    % outputs:
    % alpha = array of random angles with range [0, pi] = [0 0.5 cycle] 
    % relevant  for phase differences used in gait space calculations
    mu = mu*2*pi;                           % convert mean to radians
    % compute random array of angles in radians over range [-pi,pi]
    alpha = circ_vmrnd(mu,kappa,nsamples);  % angles in radians
    alpha(alpha<0) = 2*pi + alpha(alpha<0); % only consider positive angles
    % to limit range to desired [0, pi]
    alpha = alpha/(2*pi);                   % convert to cycles

end
%% ************************************************************************
function alpha = circ_vmrnd(theta, kappa, n)
    % alpha = circ_vmrnd(theta, kappa, n)
    %   Simulates n random angles from a von Mises distribution, with preferred 
    %   direction thetahat and concentration parameter kappa.
    %
    %   Input:
    %     [theta    preferred direction, default is 0]
    %     [kappa    width, default is 1]
    %     [n        number of samples, default is 10]
    %
    %     If n is a vector with two entries (e.g. [2 10]), the function creates
    %     a matrix output with the respective dimensionality.
    %
    %   Output:
    %     alpha     samples from von Mises distribution
    %
    %
    %   References:
    %     Statistical analysis of circular simdata, Fisher, sec. 3.3.6, p. 49
    %
    % Circular Statistics Toolbox for Matlab
    
    % By Philipp Berens and Marc J. Velasco, 2009
    % velasco@ccs.fau.edu
    
    
    % default parameter
    if nargin < 3
        n = 10;
    end
    
    if nargin < 2
        kappa = 1;
    end
    
    if nargin < 1
        theta = 0;
    end
    
    if numel(n) > 2
      error('n must be a scalar or two-entry vector!')
    elseif numel(n) == 2
      m = n;
      n = n(1) * n(2);
    end  
    
    % if kappa is small, treat as uniform distribution
    if kappa < 1e-6
        alpha = 2*pi*rand(n,1);
        return
    end
    
    % other cases
    a = 1 + sqrt((1+4*kappa.^2));
    b = (a - sqrt(2*a))/(2*kappa);
    r = (1 + b^2)/(2*b);
    
    alpha = zeros(n,1);
    for j = 1:n
      while true
          u = rand(3,1);
    
          z = cos(pi*u(1));
          f = (1+r*z)/(r+z);
          c = kappa*(r-f);
    
          if u(2) < c * (2-c) || ~(log(c)-log(u(2)) + 1 -c < 0)
             break
          end
    
          
      end
    
      alpha(j) = theta +  sign(u(3) - 0.5) * acos(f);
      alpha(j) = angle(exp(1i*alpha(j)));
    end
    
    if exist('m','var')
      alpha = reshape(alpha,m(1),m(2));
    end
end
%% ************************************************************************