function [obj,meta] = getTrialPsth(meta,obj)

edges = meta.tmin:meta.dt:meta.tmax;

% get single trial data
obj.trialpsth = zeros(numel(obj.time),numel(meta.cluNum),obj.bp.Ntrials);
obj.trialcounts = obj.trialpsth;
for i = 1:numel(meta.cluNum)
    for j = 1:obj.bp.Ntrials
        curClu = meta.cluNum(i);
        trix = j;
        
        obj.trialId(j) = trix;
        
        spkix = ismember(obj.clu{meta.probe}(curClu).trial, trix);

        N = histc(obj.clu{meta.probe}(curClu).trialtm_aligned(spkix), edges)./meta.dt;
        N = N(1:end-1);
        if size(N,2) > size(N,1)
            N = N'; % make sure N is a column vector
        end
        
        obj.trialpsth(:,i,j) = mySmooth(N,meta.smooth);
        
    end
end


end % getTrialPsth











