function obj = alignSpikes(obj,meta,params)

if strcmpi(params.alignEvent,'moveOnset')
    obj = findMoveOnset(obj); % assigns obj.bp.ev.moveOnset
end

% align spikes to params.alignEvent
for clu = 1:numel(obj.clu{meta.probe})
    event = obj.bp.ev.(params.alignEvent)(obj.clu{meta.probe}(clu).trial);
    obj.clu{meta.probe}(clu).trialtm_aligned = obj.clu{meta.probe}(clu).trialtm - event;
end

end % alignSpikes