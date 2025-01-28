%% cluster measured gait phase oscillation distances

% Description:
% takes an input file with the Nlegs-1 phase differences between adjacent
% legs, computes all pairwise distances between them in gait space and uses
% hierarchical clustering with Ward's linkage to assign them to the desired
% number of clusters

% Saves the resulting cluster assignments & computers statistics on the
% clustered data

% Input:  spreadsheet file readable by readcell with the first row = header
% and the columns in this order:
% columns 1 ... Nlegs-1 = oscillation phase differences (in cycles)
% between adjacent legs, computed as: phase i+1 - phase i
% These differences need to be computed consistently 
% for the same pairs of adjacent legs for the entire dataset.
% if you wish to use the model gait data provided, you also need to use the
% same convention for i for the models as for any measured data.
% 
% The other columns can be metadata entries that describe other aspects of
% the data, or other kinematic variables.

% Output:  spreadsheet that contains all of the input data (including
% headers) with the assigned cluster number in the first column

% References: Clustering phase difference data in gait space: Kane, Suzanne
% Amador, Brooke L. Quinn, Xuanyi Kris Wu, Sarah Y. Xi, Michael F. Ochs,
% and S. Tonia Hsieh. "Unsupervised learning reveals rapid gait adaption
% after leg loss and regrowth in spiders." bioRxiv (2025): 2025-01.
%
% gait space metric:
% Wilshin, Simon, Paul S. Shamble, Kyle J. Hovey, Ryan Harris, Andrew J.
% Spence, and S. Tonia Hsieh. "Limping following limb loss increases
% locomotor stability." Journal of Experimental Biology 221, no. 18 (2018):
% jeb174268.
%
% phase difference calculation from foot motion data: Revzen, Shai, Samuel
% A. Burden, Talia Y. Moore, Jean-Michel Mongeau, and Robert J. Full.
% "Instantaneous kinematic phase reflects neuromechanical response to
% lateral perturbations of running cockroaches." Biological cybernetics
% 107, no. 2 (2013): 179-200.

%% initialization
clear all; close all; clc;

%% get phase difference spreadsheet
[datafile, filepath] = uigetfile('*.*','Select phase difference file (1st row = header)');
cd(filepath);
% get the corefilename to use in naming other files
[~,corefilename,~] = fileparts(datafile);
% create full filename
datafilename = fullfile(filepath,datafile);

clustertree = fullfile(filepath,[corefilename,'tree.png']);

%% get filename and data from the master spreadsheet file
inputdata = readcell(datafilename);

% get number of legs
prompt = {'Number of legs in file ',corefilename,'?'};
answer = inputdlg(prompt,'',[1 40],{'4'});

% number of legs for which gaits were analyzed
Nlegs = str2num(answer{1});

if Nlegs > 8 || Nlegs < 3
    display('Only 3-8 legs supported -- see documentation for how to add option for more');
else
    
    %% get data to cluster
    
    % get the oscillation phase difference data
    % (ignore 1st row = header)
    data_to_cluster = cell2mat(inputdata(2:end,1:(Nlegs-1)));
    
    nsamples = size(data_to_cluster,1);             % number of data samples
    
    fprintf('%i samples of oscillation phase difference to cluster\n',nsamples);
    
    %% perform clustering
    
    % first compute the gait space distances between all pairs of data:
    display('computing gait space distances for hierarchical clustering')
    dist = pdist(data_to_cluster,@gaitspacedistance);
    
    % next compute Ward's linkages
    display('computing linkages for hierarchical clustering')
    hierarch_cluster_tree = linkage(dist,'ward');
    
    % show the resulting dendrogram for the linkages
    figure('Color','white','Position',[20 50 400 400]);
    hdend = dendrogram(hierarch_cluster_tree);
    
    % save the dendrogram as an image file
    saveas(gcf,clustertree);
    
    % choose number of clusters
    answer = inputdlg({'Assign data to how many clusters?'},'Number clusters',[1 40],{'3'});
    nclusters = str2num(answer{1});
    
    % cluster the data & save the cluster ID's (numbers indicating which
    % cluster each sample has been assigned to
    idxcluster = cluster(hierarch_cluster_tree,'maxclust',nclusters);
    
    % done clustering
    
    fprintf('Input data assigned to %3.0f clusters by hierarchical clustering\n',nclusters);
    fprintf('--------------------------\n');
    
    %% save input data with the cluster assignments
    
    % create output datafile
    outputfilename = fullfile(filepath,[corefilename,'_clust.csv']);
    if isfile(outputfilename)  % output spreadsheet exists -- delete
        delete outputfilename;
    end
    
    % append cluster assignments (with header) as the last column 
    outputdata = {};
    outputdata = ['cluster id';num2cell(idxcluster)];
    outputdata = [inputdata outputdata];
    
    % save input data with cluster assignments in last column
    writecell(outputdata,outputfilename);
    
    %% compute cluster circular statistics for all clusters
    display('compute and display cluster statistics');
    
    cluster_center = [];
    for i = 1:nclusters
        % get circular mean after converting to radians
        cluster_data = data_to_cluster(idxcluster == i,:);
        [clcent,clul,clll] = ...
            circ_mean(cluster_data*2*pi);
        % convert back to cycles
        cluster_center(i,:) = clcent/(2*pi);
        cluster_ul(i,:) = clul/(2*pi);
        cluster_ll(i,:) = clll/(2*pi);
        % compute the angular deviation, analogous to SD
        [AD,~] = circ_std(cluster_data*2*pi);
        ang_dev(i,:) = AD/(2*pi);  
        % compute each cluster's von Mises distribution parameters mu and kappa
        [cluster_mu(i), cluster_kappa(i)] = circ_vmpar(cluster_data*2*pi);
    end
    
    for i = 1:nclusters
        numclustersamples(i) = sum(idxcluster == i);
    end
    
    %% create a table for nice display
    cluster_statistics = array2table(cluster_center);
    
    % add column with number of samples
    cluster_statistics = addvars(cluster_statistics,numclustersamples','Before',1,'NewVariableNames','number samples');
    
    % add column with cluster number
    Cluster = [1:nclusters]';
    cluster_statistics = addvars(cluster_statistics,Cluster,'Before',1);
    
    disp(cluster_statistics);
    
    % %% print out the cluster centers
    % fprintf('cluster phase difference (cycle) mean [95%% CI] angular deviation\n');
    % for i = 1:nclusters
    %     for j = 1:size(data_to_cluster,2)
    %         fprintf('%3.3f [%3.3f, %3.3f], %3.3f ',...
    %             cluster_center(i,j),cluster_ll(i,j),cluster_ul(i,j),cluster_ul(i,j));
    %     end
    %     fprintf('\n');
    % end
    % 
    % for i = 1:nclusters
    %     fprintf('number samples in cluster %i = %5.0f\n',i,sum(idxcluster == i))
    % end
end
%% end of main program

%% ************************************************************************
function [DistAB] = gaitspacedistance(diffphaseA,diffphaseB)
% computes the gait space distances between diffphaseA (a row
% vector) and each of the Nsamples rows in diffphase B
% adapted from:
% Wilshin, S., Reeve, M. A., Haynes, G. C., Revzen, S., Koditschek, D.
% E., & Spence, A. J. (2017). Longitudinal quasi-static stability
% predicts changes in dog gait on rough terrain. Journal of
% Experimental Biology, 220(10), 1864-1874.
% Wilshin, S., Shamble, P. S., Hovey, K. J., Harris, R., Spence, A. J.,
% & Hsieh, S. T. (2018). Limping following limb loss increases
% locomotor stability. Journal of Experimental Biology, 221(18),
% jeb174268.

% VARIABLES: diffphaseA, diffphaseB correspond
% to the phi0,...phiN-1 in Eqn 1-3 in Wilshin et al., 2017, which are
% differences between each pair of phases i = 1 to N-1 (N-1 values)
% numbered as i,i-1-th pair
%
% UNITS:  We will consistently express all phases in cycles over the
% range [0,1] cycle = [0,2pi]radians.

Nlegs = size(diffphaseA,2)+1;             % number of legs used in gait
Nsamples = size(diffphaseB,1);        % number of gait samples in B

% compute the difference between diffphaseA and diffphaseB (eq. 2 & 3
% in Wilshin et al., 2018):
D_AB = diffphaseB - repmat(diffphaseA,size(diffphaseB,1),1);

% Because of periodicity, the actual distance
% along each dimension is the minimum of the difference between phases
% and 1 cycle
% Ref Wilshin et al., 2018: "We note that owing to the topology of the
% space of limb phase differences (a five-dimensional hyper-torus),
% the distances given by D are not unique, because any elements of the
% tuple which make up the arguments can be shifted by 2π (1 cycle) and
% be equivalent.  We therefore defined the distance to be the minimal
% value of D from those possible by performing such shifts."

% from Simon Wilshin's toroidalGeometry.py
D_AB = mod((D_AB+0.5),1) - 0.5;

% combine the differences along each dimension to
% get the total distance using the metric derived in Wilshin et al.,
% 2018 (eq. 6) (also see Wilshin et al., 2017)

% this metric, called gij in the paper, depends on the number of legs
% as follows:
switch Nlegs
    case 2
        DistAB = NaN;
        fprintf('gait space metric does not make sense for bipods (2 legs)\n');
        return;
    case 3              
        gij = (1/3)*[2 1;  
                     1 2];  
    case 4
        % the 2017 dog locomotion paper has an incorrect version of this 
        % tensor, which is given as g3T in Wilshin's toroidalGeometry.py code
        gij = (1/4)*[3 2 1;  
                     2 4 2;
                     1 2 3];
    case 5              
        gij = (1/5)*[4 3 2 1;
                     3 6 4 2;
                     2 4 6 3;
                     1 2 3 4];
    case 6           % g5T in Wilshin's toroidalGeometry.py code
        gij = (1/6)*[5 4 3 2 1;
                     4 8 6 4 2;
                     3 6 9 6 3;
                     2 4 6 8 4;
                     1 2 3 4 5];        
    case 7
        gij = (1/7)*[6 5  4  3  2  1;
                     5 10 8  6  4  2;
                     4 8  12 9  6  3;
                     3 6  9  12 8  4;
                     2 4  6  8  10 5];
    case 8    % g7T in Wilshin's toroidalGeometry.py code
        gij = (1/8)*[7 6  5   4   3   2   1;
                     6 12 10  8   6   4   2;
                     5 10 15  12  9   6   3;
                     4 8  12  16  12  8   4;
                     3 6  9   12  15  10  5;
                     2 4  6   8   10  12  6;
                     1 2  3   4   5   6   7];
    otherwise           % this value of Nlegs is not supported -- return to the main program & exit
        DistAB = NaN;
        fprintf('gait space metric for > 8 legs not supported -- see documentation for how to add this as an option\n',Nlegs);
        return;
end

% compute total distance using Eqn. 7 from Wilshin et al., 2018
% for each of the jth columns in D_AB corresponding to each sample in
% the input
for i = 1:Nsamples
    % compute distance squared
    dAB(i,1) = D_AB(i,:)*gij*D_AB(i,:)';
    % take the square root to get distance in cycles
    dAB(i,1) = sqrt(dAB(i,1));
end

% normalize to maximum possible distance
% maximum distances are 0.5 cycle = pi along each of Nlegs dimensions
D_AB_max = 0.5*ones(Nlegs-1,1);
dABmax = sqrt(D_AB_max'*gij*D_AB_max);
% fprintf('max distance = %f for %d legs\n',dABmax,Nlegs);
DistAB = dAB/dABmax; % fraction maximum possible dist
end
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
%     Statistical analysis of circular data_to_cluster, Fisher, sec. 3.3.6, p. 49
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
function [mu,ul,ll] = circ_mean(alpha, w, dim)
%
% mu = circ_mean(alpha, w)
%   Computes the mean direction for circular data_to_cluster.
%
%   Input:
%     alpha	sample of angles in radians
%     [w		weightings in case of binned angle data_to_cluster]
%     [dim  compute along this dimension, default is 1]
%
%     If dim argument is specified, all other optional arguments can be
%     left empty: circ_mean(alpha, [], dim)
%
%   Output:
%     mu		mean direction
%     ul    upper 95% confidence limit
%     ll    lower 95% confidence limit
%
% PHB 7/6/2008
%
% References:
%   Statistical analysis of circular data_to_cluster, N. I. Fisher
%   Topics in circular statistics, S. R. Jammalamadaka et al.
%   Biostatistical Analysis, J. H. Zar
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html

if nargin < 3
    dim = 1;
end

if nargin < 2 || isempty(w)
    % if no specific weighting has been specified
    % assume no binning has taken place
    w = ones(size(alpha));
else
    if size(w,2) ~= size(alpha,2) || size(w,1) ~= size(alpha,1)
        error('Input dimensions do not match');
    end
end

% compute weighted sum of cos and sin of angles
r = sum(w.*exp(1i*alpha),dim);

% obtain mean by
mu = angle(r);

% confidence limits if desired
if nargout > 1
    t = circ_confmean(alpha,0.05,w,[],dim);
    ul = mu + t;
    ll = mu - t;
end
end
%% ************************************************************************
function t = circ_confmean(alpha, xi, w, d, dim)
%
% t = circ_mean(alpha, xi, w, d, dim)
%   Computes the confidence limits on the mean for circular data_to_cluster.
%
%   Input:
%     alpha	sample of angles in radians
%     [xi   (1-xi)-confidence limits are computed, default 0.05]
%     [w		number of incidences in case of binned angle data_to_cluster]
%     [d    spacing of bin centers for binned data_to_cluster, if supplied
%           correction factor is used to correct for bias in
%           estimation of r, in radians (!)]
%     [dim  compute along this dimension, default is 1]
%
%   Output:
%     t     mean +- d yields upper/lower (1-xi)% confidence limit
%
% PHB 7/6/2008
%
% References:
%   Statistical analysis of circular data_to_cluster, N. I. Fisher
%   Topics in circular statistics, S. R. Jammalamadaka et al.
%   Biostatistical Analysis, J. H. Zar
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html

if nargin < 5
    dim = 1;
end

if nargin < 4 || isempty(d)
    % per default do not apply correct for binned data_to_cluster
    d = 0;
end

if nargin < 3 || isempty(w)
    % if no specific weighting has been specified
    % assume no binning has taken place
    w = ones(size(alpha));
else
    if size(w,2) ~= size(alpha,2) || size(w,1) ~= size(alpha,1)
        error('Input dimensions do not match');
    end
end

% set confidence limit size to default
if nargin < 2 || isempty(xi)
    xi = 0.05;
end

% compute ingredients for conf. lim.
r = circ_r(alpha,w,d,dim);
n = sum(w,dim);
R = n.*r;
c2 = chi2inv((1-xi),1);

% check for resultant vector length and select appropriate formula
t = zeros(size(r));

for i = 1:numel(r)
    if r(i) < .9 && r(i) > sqrt(c2/2/n(i))
        t(i) = sqrt((2*n(i)*(2*R(i)^2-n(i)*c2))/(4*n(i)-c2));  % equ. 26.24
    elseif r(i) >= .9
        t(i) = sqrt(n(i)^2-(n(i)^2-R(i)^2)*exp(c2/n(i)));      % equ. 26.25
    else
        t(i) = NaN;
        warning('Requirements for confidence levels not met.');
    end
end

% apply final transform
t = acos(t./R);

end
%% ************************************************************************
function r = circ_r(alpha, w, d, dim)
% r = circ_r(alpha, w, d)
%   Computes mean resultant vector length for circular data_to_cluster.
%
%   Input:
%     alpha	sample of angles in radians
%     [w		number of incidences in case of binned angle data_to_cluster]
%     [d    spacing of bin centers for binned data_to_cluster, if supplied
%           correction factor is used to correct for bias in
%           estimation of r, in radians (!)]
%     [dim  compute along this dimension, default is 1]
%
%     If dim argument is specified, all other optional arguments can be
%     left empty: circ_r(alpha, [], [], dim)
%
%   Output:
%     r		mean resultant length
%
% PHB 7/6/2008
%
% References:
%   Statistical analysis of circular data_to_cluster, N.I. Fisher
%   Topics in circular statistics, S.R. Jammalamadaka et al.
%   Biostatistical Analysis, J. H. Zar
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html

if nargin < 4
    dim = 1;
end

if nargin < 2 || isempty(w)
    % if no specific weighting has been specified
    % assume no binning has taken place
    w = ones(size(alpha));
else
    if size(w,2) ~= size(alpha,2) || size(w,1) ~= size(alpha,1)
        error('Input dimensions do not match');
    end
end

if nargin < 3 || isempty(d)
    % per default do not apply correct for binned data_to_cluster
    d = 0;
end

% compute weighted sum of cos and sin of angles
r = sum(w.*exp(1i*alpha),dim);

% obtain length
r = abs(r)./sum(w,dim);

% for data_to_cluster with known spacing, apply correction factor to correct for bias
% in the estimation of r (see Zar, p. 601, equ. 26.16)
if d ~= 0
    c = d/2/sin(d/2);
    r = c*r;
end

end
%% ************************************************************************
function [s,s0] = circ_std(alpha, w, d, dim)
% s = circ_std(alpha, w, d, dim)
%   Computes circular standard deviation for circular data_to_cluster
%   (equ. 26.20, Zar).
%
%   Input:
%     alpha	sample of angles in radians
%     [w		weightings in case of binned angle data_to_cluster]
%     [d    spacing of bin centers for binned data_to_cluster, if supplied
%           correction factor is used to correct for bias in
%           estimation of r]
%     [dim  compute along this dimension, default is 1]
%
%     If dim argument is specified, all other optional arguments can be
%     left empty: circ_std(alpha, [], [], dim)
%
%   Output:
%     s     angular deviation
%     s0    circular standard deviation
%
% PHB 6/7/2008
%
% References:
%   Biostatistical Analysis, J. H. Zar
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html

if nargin < 4
    dim = 1;
end

if nargin < 3 || isempty(d)
    % per default do not apply correct for binned data_to_cluster
    d = 0;
end

if nargin < 2 || isempty(w)
    % if no specific weighting has been specified
    % assume no binning has taken place
    w = ones(size(alpha));
else
    if size(w,2) ~= size(alpha,2) || size(w,1) ~= size(alpha,1)
        error('Input dimensions do not match');
    end
end

% compute mean resultant vector length
r = circ_r(alpha,w,d,dim);

s = sqrt(2*(1-r));      % 26.20
s0 = sqrt(-2*log(r));    % 26.21

end
%% ************************************************************************
function [thetahat kappa] = circ_vmpar(alpha,w,d)
    % r = circ_vmpar(alpha, w, d)
    %   Estimate the parameters of a von Mises distribution.
    %
    %   Input:
    %     alpha	sample of angles in radians
    %     [w		number of incidences in case of binned angle data]
    %     [d    spacing of bin centers for binned data, if supplied 
    %           correction factor is used to correct for bias in 
    %           estimation of r, in radians (!)]
    %
    %   Output:
    %     thetahat		preferred direction
    %     kappa       concentration parameter
    %
    % PHB 3/23/2009
    %
    % References:
    %   Statistical analysis of circular data, N.I. Fisher
    %
    % Circular Statistics Toolbox for Matlab
    
    % By Philipp Berens, 2009
    % berens@tuebingen.mpg.de
    
    alpha = alpha(:);
    if nargin < 2
      w = ones(size(alpha));
    end
    if nargin < 3
      d = 0;
    end
    
    r = circ_r(alpha,w,d);
    kappa = circ_kappa(r);
    
    thetahat = circ_mean(alpha,w);
end
%% ************************************************************************
function kappa = circ_kappa(alpha,w)
    %
    % kappa = circ_kappa(alpha,[w])
    %   Computes an approximation to the ML estimate of the concentration 
    %   parameter kappa of the von Mises distribution.
    %
    %   Input:
    %     alpha   angles in radians OR alpha is length resultant
    %     [w      number of incidences in case of binned angle data]
    %
    %   Output:
    %     kappa   estimated value of kappa
    %
    %   References:
    %     Statistical analysis of circular data, Fisher, equation p. 88
    %
    % Circular Statistics Toolbox for Matlab
    
    % By Philipp Berens, 2009
    % berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html
    
    
    alpha = alpha(:);
    
    if nargin<2
      % if no specific weighting has been specified
      % assume no binning has taken place
	    w = ones(size(alpha));
    else
      if size(w,2) > size(w,1)
        w = w';
      end 
    end
    
    N = length(alpha);
    
    if N>1
      R = circ_r(alpha,w);
    else
      R = alpha;
    end
    
    if R < 0.53
      kappa = 2*R + R^3 + 5*R^5/6;
    elseif R>=0.53 && R<0.85
      kappa = -.4 + 1.39*R + 0.43/(1-R);
    else
      kappa = 1/(R^3 - 4*R^2 + 3*R);
    end
    
    if N<15 && N>1
      if kappa < 2
        kappa = max(kappa-2*(N*kappa)^-1,0);    
      else
        kappa = (N-1)^3*kappa/(N^3+N);
      end
    end
end
%% ************************************************************************