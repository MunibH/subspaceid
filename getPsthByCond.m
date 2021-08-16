function [obj,sidmeta] = getPsthByCond(sidmeta,obj,params)

obj.edges = sidmeta.tmin:sidmeta.dt:sidmeta.tmax;
obj.time = obj.edges + sidmeta.dt/2;
obj.time = obj.time(1:end-1);

obj.condition = sidmeta.condition';

% ensure spikes are aligned to go cue
for clu = 1:numel(obj.clu{sidmeta.probe})
    event = obj.bp.ev.(params.evName)(obj.clu{sidmeta.probe}(clu).trial);
    obj.clu{sidmeta.probe}(clu).trialtm_aligned = obj.clu{sidmeta.probe}(clu).trialtm - event;
end


% get psths by condition
obj.psth = zeros(numel(obj.time),numel(sidmeta.cluNum),numel(obj.condition));
for i = 1:numel(sidmeta.cluNum)
    for j = 1:numel(obj.condition)
        curClu = sidmeta.cluNum(i);
        trix = sidmeta.trialNum{j};
        spkix = ismember(obj.clu{sidmeta.probe}(curClu).trial, trix);

        N = histc(obj.clu{sidmeta.probe}(curClu).trialtm_aligned(spkix), obj.edges);
        N = N(1:end-1);

        obj.psth(:,i,j) = mySmooth(N./numel(trix)./sidmeta.dt, sidmeta.smooth);  % trial-averaged separated by trial type
    end
end

% get single trial data
obj.trialpsth = zeros(numel(obj.time),numel(sidmeta.cluNum),obj.bp.Ntrials);
for i = 1:numel(sidmeta.cluNum)
    for j = 1:obj.bp.Ntrials
        curClu = sidmeta.cluNum(i);
        trix = j;
        spkix = ismember(obj.clu{sidmeta.probe}(curClu).trial, trix);

        N = histc(obj.clu{sidmeta.probe}(curClu).trialtm_aligned(spkix), obj.edges);
        N = N(1:end-1);
        if size(N,2) > size(N,1)
            N = N'; % make sure N is a column vector
        end
        
        obj.trialpsth(:,i,j) = mySmooth(N./numel(trix)./sidmeta.dt, sidmeta.smooth);

    end
end

end % getPsthByCond











