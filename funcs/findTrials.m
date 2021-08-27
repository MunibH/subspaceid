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