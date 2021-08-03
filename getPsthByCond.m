function [obj,meta] = getPsthByCond(meta,obj)

obj.edges = meta.tmin:meta.dt:meta.tmax;
obj.time = obj.edges + meta.dt/2;
obj.time = obj.time(1:end-1);

obj.condition = meta.condition';

% get list of clusters
cluQuality = {obj.clu{meta.probe}(:).quality}';
meta.cluNum = findClusters(cluQuality, meta.quality);

% get list of trials
meta.trialNum = findTrials(obj, obj.condition);

% ensure spikes are aligned to go cue
for clu = 1:numel(obj.clu{meta.probe})
    event = obj.bp.ev.goCue(obj.clu{meta.probe}(clu).trial);
    obj.clu{meta.probe}(clu).trialtm_aligned = obj.clu{meta.probe}(clu).trialtm - event;
end


% get psths by condition
obj.psth = zeros(numel(obj.time),numel(meta.cluNum),numel(obj.condition));
for i = 1:numel(meta.cluNum)
    for j = 1:numel(obj.condition)
        curClu = meta.cluNum(i);
        trix = meta.trialNum{j};
        spkix = ismember(obj.clu{meta.probe}(curClu).trial, trix);

        N = histc(obj.clu{meta.probe}(curClu).trialtm_aligned(spkix), obj.edges);
        N = N(1:end-1);

        obj.psth(:,i,j) = mySmooth(N./numel(trix)./meta.dt, meta.smooth);  % trial-averaged separated by trial type
    end
end

end % getPsthByCond

%%

function idx = findClusters(qualityList, qualities)
% find idx where qualityList contains at least one of the patterns in
% qualities

% handle unlabeled cluster qualities
for i = 1:numel(qualityList)
    if isempty(qualityList{i})
        qualityList(i) = {'nan'};
    end
end

[~,mask] = patternMatchCellArray(qualityList, qualities, 'any');

idx = find(mask);
end % findClusters

function trialNums = findTrials(obj, conditions)

varnames = getStructVarNames(obj);
for i = 1:numel(varnames)
    eval([varnames{i} ' = obj.bp.' varnames{i} ';']);
    
    if eval(['numel(' varnames{i} ')==obj.bp.Ntrials && isrow(' varnames{i} ')'])
        eval([varnames{i} '=' varnames{i} ''';']);
    end
end

mask = zeros(obj.bp.Ntrials, numel(conditions));
for i = 1:numel(conditions)
    mask(:,i) = eval(conditions{i}); 
    trialNums{i} = find(mask(:,i));
end

end % findTrials







