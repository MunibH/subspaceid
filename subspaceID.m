clear,clc,close all

if ispc
    pth = 'C:\Code\subspace-id';
elseif ismac
    pth = '/Users/Munib/Documents/Economo-Lab/code/subspaceid';
end
addAllPaths(pth);

%% TODO
%  fix bug in optimization when you change literally anything
%  var explained - opt method
%  move onset
%  handle multiple probes of data
%  handle more than 2 conditions
%  kaufman method (regression)
%  dpca
%  psid
%  gpfa for single trials

%% SET RUN PARAMS

params.doPSTH              = false; % get psths by condition and save meta data from above
params.evName              = 'goCue';
params.sav                 = 0;    % save obj with psths (just need to do this once)

params.method.optimization = true;   % elsayed method
params.method.maxdiff      = false;  % new method mike and chand came up with
params.method.regression   = false;  % kaufman method

params.conditions          = [1 , 3]; % which conditions to use in analysis (only 2 rn)

params.varToExplain        = 90;    % sets dimensionality of null and potent space

% for 3D plot
params.dims.potent         = [1,2]; % potent dims to plot by default
params.dims.null           = [1];   % null dims to plot by default

%% SET METADATA
% experiment meta data
meta.datapth = fullfile(pth,'data');
meta.anm = 'JEB7';
meta.date = '2021-04-29';
meta.probe = 1;

% analysis meta data
meta.tmin = -2.2; % (s) relative to go cue
meta.tmax = 3;  % (s) relative to go cue
meta.dt = 0.005;
meta.smooth = 200; % for params.doPSTH
meta.prepEpoch = [-2.2, -0.15]; % (s) relative to go cue
meta.moveEpoch = [0.15, 1.15]; % (s) relative to go cue

% conditions (i.e. trials to look at)
meta.condition(1) = {'R&hit&~stim.enable&autowater.nums==2'}; % right hits, no stim, aw off
meta.condition(2) = {'L&hit&~stim.enable&autowater.nums==2'}; % left hits, no stim, aw off
meta.condition(3) = {'R&hit&~stim.enable&autowater.nums==1'}; % right hits, no stim, aw on
meta.condition(4) = {'L&hit&~stim.enable&autowater.nums==1'}; % left hits, no stim, aw on
% meta.condition(5) = {'hit&~stim.enable&autowater.nums==2'};   % hits, no stim, aw off
% meta.condition(6) = {'hit&~stim.enable&autowater.nums==1'};   % hits, no stim, aw on

% clusters (these qualities are included)
meta.quality = {'Fair','Good','Great','Excellent','single','multi'}; 

%% SET RUN PARAMS

params.doPSTH              = false; % get psths by condition and save meta data from above

params.method.optimization = true;  % elsayed method
params.method.maxdiff      = false;  % new method mike and chand came up with
params.method.regression   = false; % kaufman method

params.conditions          = [1,2]; % which conditions to use in analysis (only 2 rn)

params.varToExplain        = 90;    % sets dimensionality of null and potent space

% for 3D plot
params.dims.potent         = [1,2]; % potent dims to plot by default
params.dims.null           = [1];   % null dims to plot by default


%% LOAD DATA
[obj,meta] = loadDataObj(meta);

%% get trials and clusters to use

% get list of clusters
cluQuality = {obj.clu{meta.probe}(:).quality}';
meta.cluNum = findClusters(cluQuality, meta.quality);

% get list of trials
meta.trialNum = findTrials(obj, obj.condition);

%% PSTHs
if params.doPSTH
    [obj, meta] = getPsthByCond(meta,obj,params);
    if params.sav
        save(fullfile(meta.datapth, meta.datafn), 'obj', '-v7.3');
    end
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
addpath(genpath(fullfile(pth,'preprocess')))
addpath(genpath(fullfile(pth,'optimization')))
addpath(genpath(fullfile(pth,'regression')))
addpath(genpath(fullfile(pth,'maxdiff_pca')))

end % addAllPaths

function [obj,meta] = loadDataObj(meta)
contents = dir(meta.datapth);
contents = {contents.name}';

strToFind = {'data_structure' , meta.anm, meta.date};

[meta.datafn,~] = patternMatchCellArray(contents, strToFind, 'all');

dat = load(fullfile(meta.datapth, meta.datafn));
obj = dat.obj;

end % loadRawDataObj

function rez = subspaceIDWithMethod(meta, obj, method, params)
switch method
    case 'optimization'
        rez = subspaceid_optimization(meta, obj, params);
    case 'regression'
        rez = subspaceid_regression(meta, obj, params);
    case 'maxdiff'
        rez = subspaceid_maxdiff(meta, obj, params);
end
end % subspaceIDWithMethod












