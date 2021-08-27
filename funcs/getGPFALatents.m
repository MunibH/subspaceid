function latents = getGPFALatents(obj,meta,params)

dat = formatDat(obj);

% GPFA
resultsPth = 'C:\Code\subspace-id\gpfa\mat_results\';
runIdx = getNextRunIdx(resultsPth); 
method = 'gpfa';
xDim = 4; 
kernSD = 30; % 30 (default) - for two-stage methods
segLength = 20; % 20 (default) - didn't notice differences in changing this
binWidth = 20; % 20 (default) - differences in changing this param depends on dt
result = neuralTraj(runIdx, dat(1:40), 'datFormat', 'seq',...
                    'method', method, 'xDim', xDim,... 
                    'kernSDList', kernSD, 'segLength', segLength,...
                    'binWidth', binWidth);

[~, seqTrain] = postprocess(result);

latents = nan(seqTrain(1).T,numel(seqTrain),xDim); % (time,trials,dims)
for n = 1:numel(seqTrain) % for every trial
    for k = 1:xDim % for every latent dim
        latents(:,n,k) = seqTrain(n).xorth(k,:);
    end
end

% MAKE SOME BETTER PLOTTING FUNCTIONS (A GUI LIKE THE OTHERS WOULD BE COOL)

redTrials = meta.trialNum{1};

nPlotMax = 100;
plotEachDimVsTime(seqTrain, 'xorth', obj.time,'nPlotMax',nPlotMax,...
                  'redTrials',redTrials);
nPlotMax3D = 20;
plot3D(seqTrain, 'xorth', 'dimsToPlot', [1 2 4], 'nPlotMax',nPlotMax3D,...
                  'redTrials',redTrials);


end % getGPFALatents


function dat = formatDat(obj)
% format data for gpfa
% if datFormat is 'seq', nth entry has fields
%                       trialId      -- unique trial identifier  
%                       T (1 x 1)    -- number of timesteps
%                       y (yDim x T) -- continuous valued data 
%                                       (Eg: binned spike counts)
for i = 1:size(obj.trialcounts,3)
    dat(i).trialId = obj.trialId(i);
    dat(i).T       = numel(obj.time);
    dat(i).y = obj.trialcounts(:,:,i)';
end
end % formatDat

function runIdx = getNextRunIdx(resultsPth)
    contents = dir(resultsPth);
    contents = contents([contents.isdir]');
    dirs = {contents.name}';
    
    runIdx = 0;
    for i = 1:numel(dirs)
        if ~contains(dirs{i},'run')
            continue
        end
        curRunIdx = str2double(dirs{i}(4:end));
        if curRunIdx > runIdx
            runIdx = curRunIdx;
        end
    end
    runIdx = runIdx + 1;
end % getNextRunIdx








