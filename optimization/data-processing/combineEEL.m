% script to combine psths from multiple animals and/or sessions - EEL

data_pth = '/Users/Munib/Documents/Economo-Lab/subspace-id/data/2afc/ALM/';
files = dir([data_pth '*.mat']);

combined_psth = cell(1,numel(files));
for file = 1:numel(files)
    load(fullfile(data_pth,files(file).name))
    
    trials = getTrialsToUse(obj);
    clusters = getClustersToUse(obj);
    
    dt = 0.005; % 50 ms, same as mike's data
    time = (mode(obj.bp.ev.delay)-3):dt:(mode(obj.bp.ev.goCue)+3);
    
    psth = get_psth(obj,time,trials,clusters);
    combined_psth{file} = psth;
    
end

tag.psth = [combined_psth{1};combined_psth{2}];
tag.time = time - mode(obj.bp.ev.goCue);

save('alldat.mat','tag');

function trialsToUse = getTrialsToUse(obj)
% narrow down the trials to use to trials where delay and go cue occur at
% most frequent time and where there were successful licks. 
    % get trials where delay epoch starts at mode of obj.bp.ev.delay
    % and go cue starts at mode of obj.bp.ev.goCue
    delay = obj.bp.ev.delay; % (nTrials,1) double 
    goCue = obj.bp.ev.goCue; % (nTrials,1) double 
    modeDelay = mode(delay);
    modeGoCue = mode(goCue);
    trialsMask = (delay==modeDelay) & (goCue==modeGoCue);    
    % get successful lick trials
    right_hit_trials = obj.bp.hit & obj.bp.R;
    left_hit_trials = obj.bp.hit & obj.bp.L;
    trialsMask_right = trialsMask & right_hit_trials;
    trialsMask_left = trialsMask & left_hit_trials;    
    % get trials to use
    trials_right = find(trialsMask_right);
    trials_left = find(trialsMask_left);        
    % {right,left}
    trialsToUse = {trials_left, trials_right};
end

function clustersToUse = getClustersToUse(obj)
% get clusters of quality not equal to 'Poor' to use for analysis
    numClusters = numel(obj.clu{2});
    clustersToUse = zeros(1,numClusters);
    ct = 1;
    for cluster = 1:numClusters
        if ~strcmp(obj.clu{1}(cluster).quality,'Poor')
            clustersToUse(ct) = cluster;
            ct = ct + 1;
        end
    end
    % remove unused space
    clustersToUse(clustersToUse==0) = [];    
end

function psth = get_psth(obj,time,trials,clusters)
    numClusters = length(clusters);
    psth = cell(numClusters, 1); % {left, right}
    for cluIdx = 1:numClusters
        psth{cluIdx} = cell(1,numel(trials)); % {left,right})
        data = obj.clu{1}(clusters(cluIdx));
        for tt = 1:numel(trials)
            psth{cluIdx}{tt} = zeros(length(time),numel(trials{tt}));
            for trialIdx = 1:numel(trials{tt})
                trialtm = data.trialtm(ismember(data.trial, trials{tt}(trialIdx)));
                psth{cluIdx}{tt}(:,trialIdx) = histc(trialtm,time);
            end
        end
    end
end
