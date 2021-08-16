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