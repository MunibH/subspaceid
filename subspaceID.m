clear,clc,close all

if ispc
    pth = 'C:\Code\subspace-id';
elseif ismac
    pth = '/Users/Munib/Documents/Economo-Lab/code/subspaceid';
end

% addAllPaths(pth);
addpath(genpath(pwd))

%% TODO
%  form null and potent spaces out of movement and nonmovement indices
%  handle multiple probes of data
%  handle more than 2 conditions
%  plotting functions for gpfa latents

% OTHER METHODS:
%  kaufman method (regression)
%  dpca
%  psid
%  gpfa for single trials
%  bayesian gpfa for single trials and allows inference
%  lfads for single trials and fitting a dynamical system (predictive)

%% SET RUN PARAMS

params.method.optimization = true;   % elsayed method
params.method.maxdiff      = false;  % new method mike and chand came up with
params.method.regression   = false;  % kaufman method
params.method.psid         = false;

params.singleTrialAnalysis = true;
params.lowFR               = 5; % when doing single trial analysis, remove clusters with avg firing rate < params.lowFR

params.alignEvent          = 'moveOnset'; %'goCue'

params.conditions          = [1 , 2]; % which conditions to use in analysis (only 2 rn)

params.varToExplain        = 75;    % sets dimensionality of null and potent space

% for 3D plot
params.dims.potent         = [1,2]; % potent dims to plot by default
params.dims.null           = [1];   % null dims to plot by default

params.cols = {[98, 189, 65], [255, 57, 90]};
params.cols = cellfun(@(x) 1/255.*x, params.cols, 'UniformOutput',false);

%% SET METADATA
% experiment meta data
meta.datapth = fullfile(pth,'data');
meta.anm = 'JEB7';
meta.date = '2021-04-29';
meta.datafn = findDataFn(meta);
meta.probe = 1;

% analysis meta data
meta.tmin = -2.2; % (s) relative to params.evName
meta.tmax = 3;  % (s) relative to params.evName
meta.dt = 0.005;

meta.smooth = 200; % for params.doPSTH
meta.prepEpoch = [-1.5, -0.15]; % (s) relative to params.evName
meta.moveEpoch = [0.15, 1.15]; % (s) relative to params.evName

meta.smooth = 100; % for params.doPSTH
meta.prepEpoch = [-2.2, -0.15]; % (s) relative to params.evName
meta.moveEpoch = [0.15, 1.15]; % (s) relative to params.evName

% conditions (i.e. trials to look at)
meta.condition(1) = {'R&hit&~stim.enable&autowater.nums==2'}; % right hits, no stim, aw off
meta.condition(2) = {'L&hit&~stim.enable&autowater.nums==2'}; % left hits, no stim, aw off
meta.condition(3) = {'R&hit&~stim.enable&autowater.nums==1'}; % right hits, no stim, aw on
meta.condition(4) = {'L&hit&~stim.enable&autowater.nums==1'}; % left hits, no stim, aw on
% meta.condition(5) = {'hit&~stim.enable&autowater.nums==2'};   % hits, no stim, aw off
% meta.condition(6) = {'hit&~stim.enable&autowater.nums==1'};   % hits, no stim, aw on

% clusters (these qualities are included)
meta.quality = {'Fair','Good','Great','Excellent','single','multi'}; 


%% LOAD DATA
dat = load(fullfile(meta.datapth, meta.datafn));
obj = dat.obj;

%% get trials and clusters to use

% get list of clusters
cluQuality = {obj.clu{meta.probe}(:).quality}';
meta.cluNum = findClusters(cluQuality, meta.quality);

% get list of trials
obj.condition = meta.condition(params.conditions)';
meta.trialNum = findTrials(obj, obj.condition);

%% PSTHs and GPFA Latents 
[obj, meta] = getPsthByCond(meta,obj,params);

%% GPFA Latents for single trial analysis
if params.singleTrialAnalysis
    [obj,meta] = removeLowFRClusters(obj,meta,params);
    obj.gpfalatents = getGPFALatents(obj,meta,params);
end

%% ANALYSIS METHODS
methods = fieldnames(params.method);
ct = 1;
for i = 1:numel(methods)
    if params.method.(methods{i})
        rez{ct} = subspaceIDWithMethod(meta, obj, methods{i}, params);
        rez{ct}.method = methods{i};
        ct = ct + 1;
    end
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function addAllPaths(pth)
addpath(genpath(fullfile(pth,'utils')))
addpath(genpath(fullfile(pth,'functions')))
addpath(genpath(fullfile(pth,'plotting')))
addpath(genpath(fullfile(pth,'optimization')))
addpath(genpath(fullfile(pth,'regression')))
addpath(genpath(fullfile(pth,'maxdiff_pca')))
addpath(genpath(fullfile(pth,'psid')))

end % addAllPaths

function fn = findDataFn(meta)
contents = dir(meta.datapth);
contents = {contents.name}';

strToFind = {'data_structure' , meta.anm, meta.date};

[fn,~] = patternMatchCellArray(contents, strToFind, 'all');

end % loadRawDataObj

function rez = subspaceIDWithMethod(meta, obj, method, params)
switch method
    case 'optimization'
        rez = subspaceid_optimization(meta, obj, params);
    case 'regression'
        rez = subspaceid_regression(meta, obj, params);
    case 'maxdiff'
        rez = subspaceid_maxdiff(meta, obj, params);
    case 'psid'
        rez = subspaceid_psid(meta, obj, params);
end
end % subspaceIDWithMethod












