clear,clc,close all

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

params.method.optimization               = false;   % elsayed method
params.method.optimization_moveLabeled   = false;   % elsayed method, but null and potent are restricted to temporal epochs (TODO)
params.method.maxdiff                    = true;   % new method
% params.method.regression     = false;   % kaufman method (not working yet)
% params.method.psid           = false;   % pref subspace identification (not working yet)


params.xDim                = 6;       % number of latent dims
params.varToExplain        = 75;    % sets dimensionality of null and potent space

params.alignEvent          = 'goCue'; % 'goCue' or 'moveOnset'

% conditions/trial types to PSTHs by
params.condition(1)         = {'R&hit&~autowater'};    % right hits, no stim, aw off
params.condition(end+1)     = {'L&hit&~autowater'};    % left hits, no stim, aw off
params.condition(end+1)     = {'R&hit&autowater'};    
params.condition(end+1)     = {'L&hit&autowater'};    
params.condition(end+1)     = {'hit|miss|no'};    


params.condToUse = [1 2];

% for 3D plot
params.dims.potent         = [1,2]; % potent dims to plot by default
params.dims.null           = [1];   % null dims to plot by default

params.lowFR               = 0;      % when doing single trial analysis, remove clusters with avg firing rate < params.lowFR

params.cols = getColors();

params.probe = 2;
params.probeArea = 'ALM';

params.tmin = -2.5;
params.tmax = 1.5;
params.dt = 1/400;

params.smooth = 41;

params.quality = {'all'};

% for elsayed optimization method
% define time regions in trial that are considered prep and move times
params.prepEpoch = [-(abs(params.tmin)), -0.05]; % (s) relative to alignEvent
params.moveEpoch = [0.05, abs(params.tmax)]; % (s) relative to alignEvent

%% LOAD DATA

meta.datapth = '/Volumes/MUNIB_SSD/Experiments';
meta.datafn  = 'data_structure_EKH3_2021-08-11';

[meta,params,obj] = loadAndProcessData(meta,params);

%% label move or non-move
[obj.movix,obj.movtime] = getMoveIdx(obj,params);

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


function rez = subspaceIDWithMethod(meta, obj, method, params)
switch method
    case 'optimization'
        rez = subspaceid_optimization(obj, params);
    case 'optimization_moveLabeled'
        rez = subspaceid_optimization_moveLabeled(obj,params);
    case 'regression'
        rez = subspaceid_regression(meta, obj, params);
    case 'maxdiff'
            rez = subspaceid_maxdiff(meta, obj, params);
    case 'psid'
        rez = subspaceid_psid(meta, obj, params);
end
end % subspaceIDWithMethod












