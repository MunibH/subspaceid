function [meta,params,obj] = loadAndProcessData(meta,params)

%% LOAD DATA

dat = load(fullfile(meta.datapth,meta.datafn));
obj = dat.obj;

%% STANDARD ROUTINES

% find trials to use
params.trialid = findTrials(obj, params.condition);

% find clusters to use
params.cluid = findClusters({obj.clu{params.probe}(:).quality}', params.quality);

% align spikes in every cluster to an event
obj = alignSpikes(obj,params);

% get trial avg psth and single trial data
obj = getSeq(obj,params);

% remove low fr clusters
[obj, params] = removeLowFRClusters(obj,params);

% get mean fr for each cluster/trial type during presample (baseline fr)
obj.presampleFR = baselineFR(obj,params);


end 
