function obj = getSeq(obj,params)

edges = params.tmin:params.dt:params.tmax;
obj.time = edges + params.dt/2;
obj.time = obj.time(1:end-1);

% get psths by condition
obj.psth = zeros(numel(obj.time),numel(params.cluid),numel(params.condition));
for i = 1:numel(params.cluid)
    curClu = params.cluid(i);
    for j = 1:numel(params.condition)
        trix = params.trialid{j};
        spkix = ismember(obj.clu{params.probe}(curClu).trial, trix);

        N = histc(obj.clu{params.probe}(curClu).trialtm_aligned(spkix), edges);
        N = N(1:end-1);

        obj.psth(:,i,j) = mySmooth(N./numel(trix)./params.dt, params.smooth);  % trial-averaged separated by trial type
    end
end

% get single trial data
obj.trialpsth = zeros(numel(obj.time),numel(params.cluid),obj.bp.Ntrials);
for i = 1:numel(params.cluid)
    curClu = params.cluid(i);
    for j = 1:obj.bp.Ntrials
                
        spkix = ismember(obj.clu{params.probe}(curClu).trial, j);

        N = histc(obj.clu{params.probe}(curClu).trialtm_aligned(spkix), edges);
        N = N(1:end-1);
        if size(N,2) > size(N,1)
            N = N'; % make sure N is a column vector
        end
        
        obj.trialpsth(:,i,j) = mySmooth(N./params.dt,params.smooth);

    end
end


end % getSeq