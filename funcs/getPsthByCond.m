function [obj,meta] = getPsthByCond(meta,obj,params)

edges = meta.tmin:meta.dt:meta.tmax;
obj.time = edges + meta.dt/2;
obj.time = obj.time(1:end-1);

if strcmpi(params.alignEvent,'moveOnset')
    obj = findMoveOnset(obj); % assigns obj.bp.ev.moveOnset
end

% align spikes to params.alignEvent
for clu = 1:numel(obj.clu{meta.probe})
    event = obj.bp.ev.(params.alignEvent)(obj.clu{meta.probe}(clu).trial);
    obj.clu{meta.probe}(clu).trialtm_aligned = obj.clu{meta.probe}(clu).trialtm - event;
end

if contains(meta.datafn,'JEB7')
    tr = find(obj.bp.hit&~obj.bp.early&(obj.bp.autowater.nums'==2)); % JEB7
else
    tr = find(obj.bp.hit&~obj.bp.early&~obj.bp.autowater); % EKH2
end
Ntrials = numel(tr);

% get psths by condition
obj.psth = zeros(numel(obj.time),numel(meta.cluNum),numel(obj.condition));
for i = 1:numel(meta.cluNum)
    curClu = meta.cluNum(i);
    for j = 1:numel(obj.condition)
        trix = meta.trialNum{j};
        spkix = ismember(obj.clu{meta.probe}(curClu).trial, trix);

        N = histc(obj.clu{meta.probe}(curClu).trialtm_aligned(spkix), edges);
        N = N(1:end-1);

        obj.psth(:,i,j) = mySmooth(N./numel(trix)./meta.dt, meta.smooth);  % trial-averaged separated by trial type
    end
end

% get single trial data
nCondTrials = sum(cellfun(@numel, meta.trialNum));
obj.trialpsth = zeros(numel(obj.time),numel(meta.cluNum),nCondTrials);
obj.trialcounts = obj.trialpsth;
for i = 1:numel(meta.cluNum)
    trialct = 1;
    for j = 1:obj.bp.Ntrials
        curClu = meta.cluNum(i);
        trix = j;
        if ~(ismember(trix,meta.trialNum{1}) || ismember(trix,meta.trialNum{2}))
            continue
        end
        
        obj.trialId(trialct) = trix;
        
        spkix = ismember(obj.clu{meta.probe}(curClu).trial, trix);

        N = histc(obj.clu{meta.probe}(curClu).trialtm_aligned(spkix), edges);
        N = N(1:end-1);
        if size(N,2) > size(N,1)
            N = N'; % make sure N is a column vector
        end
        
        obj.trialcounts(:,i,trialct) = N;
        obj.trialpsth(:,i,trialct) = mySmooth(N./numel(trix)./meta.dt, meta.smooth);
        
        trialct = trialct + 1;

    end
end

end % getPsthByCond











