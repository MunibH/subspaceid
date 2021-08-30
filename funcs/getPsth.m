function [obj,meta] = getPsth(meta,obj)

edges = meta.tmin:meta.dt:meta.tmax;
obj.time = edges + meta.dt/2;
obj.time = obj.time(1:end-1);

% get psths over all trials
obj.allpsth = zeros(numel(obj.time),numel(meta.cluNum));
for i = 1:numel(meta.cluNum)
    curClu = meta.cluNum(i);
    trix = 1:obj.bp.Ntrials;
    spkix = ismember(obj.clu{meta.probe}(curClu).trial, trix);
    
    N = histc(obj.clu{meta.probe}(curClu).trialtm_aligned(spkix), edges);
    N = N(1:end-1);
    
    obj.allpsth(:,i) = mySmooth(N./numel(trix)./meta.dt, meta.smooth);  % trial-averaged over all trials
end

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



end % getPsth











