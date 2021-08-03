function subspace_id()

addAllPaths();

%% SET METADATA
% experiment meta data
meta.datapth = 'C:\Code\subspace-id\data';
meta.anm = 'JEB7';
meta.date = '2021-04-29';
meta.probe = 1;

% analysis meta data
meta.tmin = -4; % (s)relative to go cue
meta.tmax = 4;  % (s) relative to go cue
meta.dt = 0.005;
meta.smooth = 200;
meta.prepEpoch = [-1.15, -0.15]; % (s) relative to go cue
meta.moveEpoch = [0.15, 1.15]; % (s) relative to go cue

% conditions (i.e. trials to look at)
meta.condition(1) = {'R&hit&~stim.enable&autowater.nums==2'}; % right hits, no stim, aw off
meta.condition(2) = {'L&hit&~stim.enable&autowater.nums==2'}; % left hits, no stim, aw off
meta.condition(3) = {'R&hit&~stim.enable&autowater.nums==1'}; % right hits, no stim, aw on
meta.condition(4) = {'L&hit&~stim.enable&autowater.nums==1'}; % left hits, no stim, aw on

% clusters (these qualities are included)
meta.quality = {'Fair','Good','Great','Excellent','single','multi'}; 

%% SET RUN PARAMS
params.formatPSTH          = true;
params.method.optimization = false; % elsayed method
params.method.regression   = false; % kaufman method
params.method.maxdiff      = true; % new method mike and chand came up with

params.conditions          = [1,2]; % which conditions to use in analysis (only 2 rn)

%% LOAD DATA
[obj,meta] = loadRawDataObj(meta);

%% FORMAT PSTHs
if params.formatPSTH
    [obj, meta] = getPsthByCond(meta,obj);
    save(fullfile(meta.datapth, meta.datafn), 'obj', 'meta', '-v7.3');
end

%% ANALYSIS METHODS
methods = fieldnames(params.method);
for i = 1:numel(methods)
    if params.method.(methods{i})
        subspaceIDWithMethod(meta, obj, methods{i}, params);
    end
end

end % subspace_id

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function addAllPaths()
addpath(genpath('C:\Code\subspace-id\utils'))
addpath(genpath('C:\Code\subspace-id\preprocess'))
addpath(genpath('C:\Code\subspace-id\optimization'))
addpath(genpath('C:\Code\subspace-id\regression'))
addpath(genpath('C:\Code\subspace-id\maxdiff_pca'))
end % addAllPaths

function [obj, meta] = loadRawDataObj(meta)
contents = dir(meta.datapth);
contents = {contents.name}';

strToFind = {'data_structure' , meta.anm, meta.date};

[meta.datafn,~] = patternMatchCellArray(contents, strToFind, 'all');

dat = load(fullfile(meta.datapth, meta.datafn));
obj = dat.obj;
if isfield(dat,'meta')
    meta = dat.meta;
end

end % loadRawDataObj

function subspaceIDWithMethod(meta, obj, method, params)
switch method
    case 'optimization'
        doOptim(meta, obj, params)
    case 'regression'
        doRegress(meta, obj, params)
    case 'maxdiff'
        subspaceid_maxdiff(meta, obj, params)
end
end % subspaceIDWithMethod












