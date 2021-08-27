function [obj,meta] = removeLowFRClusters(obj,meta,params)
% % Remove low-firing rate units, e.g., all those firing less than 5
%   spikes per second on average across all trials.
%   
%   The fitted observation noise (diagonal element of R) for a
%   low-firing rate unit will be small, so the neural trajectory may
%   show a deflection each time the neuron spikes.

meanFRs = mean(mean(obj.psth,3));
use = meanFRs > params.lowFR;

% remove low fr clusters
meta.cluNum = meta.cluNum(use);
obj.psth = obj.psth(:,use,:);
obj.trialpsth = obj.trialpsth(:,use,:);
obj.trialcounts = obj.trialcounts(:,use,:);

end % removeLowFRClusters






