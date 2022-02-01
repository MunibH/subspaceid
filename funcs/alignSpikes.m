function obj = alignSpikes(obj,params)

if strcmp(params.alignEvent,'moveOnset')
    obj.bp.ev.(params.alignEvent) = findMoveOnset(obj);
end

if strcmp(params.alignEvent,'firstLick')
    % get first lick time for left and right licks
    temp = obj.bp.ev.lickL;
    idx = ~cellfun('isempty',temp);
    outL = zeros(size(temp));
    outL(idx) = cellfun(@(v)v(1),temp(idx));
    temp = obj.bp.ev.lickR;
    idx = ~cellfun('isempty',temp);
    outR = zeros(size(temp));
    outR(idx) = cellfun(@(v)v(1),temp(idx));
    firstLick = zeros(size(temp));
    % firstLick = min(outL,outR), except when outL||outR == 0
    outL(outL==0) = nan;
    outR(outR==0) = nan;
    firstLick = nanmin(outL,outR);
    firstLick(isnan(firstLick)) = 0;
    obj.bp.ev.(params.alignEvent) = firstLick;
end
if strcmp(params.alignEvent,'lastLick')
    % get last lick time for left and right licks
    temp = obj.bp.ev.lickL;
    idx = ~cellfun('isempty',temp);
    outL = zeros(size(temp));
    outL(idx) = cellfun(@(v)v(end),temp(idx));
    temp = obj.bp.ev.lickR;
    idx = ~cellfun('isempty',temp);
    outR = zeros(size(temp));
    outR(idx) = cellfun(@(v)v(end),temp(idx));
    lastLick = zeros(size(temp));
    % lastLick = max(outL,outR), except when outL||outR == 0
    outL(outL==0) = nan;
    outR(outR==0) = nan;
    lastLick = nanmax(outL,outR);
    lastLick(isnan(lastLick)) = 0;
    obj.bp.ev.(params.alignEvent) = lastLick;
end

% align spikes to params.alignEvent
for clu = 1:numel(obj.clu{params.probe})
    event = obj.bp.ev.(params.alignEvent)(obj.clu{params.probe}(clu).trial);
    obj.clu{params.probe}(clu).trialtm_aligned = obj.clu{params.probe}(clu).trialtm - event;
end

end % alignSpikes
