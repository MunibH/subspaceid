function rez = subspaceid_maxdiff(meta, obj, params)

%% PREPROCESS
obj.psth = obj.psth(:,:,params.conditions);
rez.time = obj.time;
rez.psth = obj.psth;
rez.trialpsth = obj.trialpsth;
rez.condition = meta.condition(params.conditions);

[obj, meta] = preprocess_maxdiff(meta, obj);

%% PREP and MOVE EPOCHS
% prepix and moveix corresponds to time idx for each epoch
[rez.prepix, rez.moveix] = getEpochs(obj.time, meta.prepEpoch, meta.moveEpoch);

%% METHOD 1 ( N*(I-WW') )

%% FIND NULL MODES
psth = obj.psth;

preppsth = psth(rez.prepix,:,:);
preppsth = [preppsth(:,:,1) ; preppsth(:,:,2)]; % (ct,n)

% pca
[pcs,explained] = myPCA(preppsth);

% find number of pcs to explain 95% of variance
rez.dPrep = numComponentsToExplainVariance(explained, params.varToExplain);
rez.Qnull = pcs(:,1:rez.dPrep);
rez.Qnull_ve = explained(1:rez.dPrep);

%% PROJECT OUT PROMINENT NULL MODES
modesToKeep = eye(size(pcs,1)) - (rez.Qnull*rez.Qnull');

proj = nan(size(psth));
for i = 1:size(psth,3)
    proj(:,:,i) = psth(:,:,i) * modesToKeep;
end

%% FIND POTENT MODES AS PCs OF MOVE EPOCH
% moveproj = proj(rez.moveix,:,:); % only use move epoch for potent mode ID
moveproj = proj; % use all leftover data for potent mode ID (seems more right)
moveproj = [moveproj(:,:,1) ; moveproj(:,:,2)]; % (ct,n)

% pca
[pcs,explained] = myPCA(moveproj);

% find number of pcs to explain 95% of variance
rez.dMove = numComponentsToExplainVariance(explained, params.varToExplain);
rez.Qpotent = pcs(:,1:rez.dMove);
rez.Qpotent_ve = explained(1:rez.dMove);

%% PLOTS
cols = {[0,0,1],[1,0,0]};
plotLatents(obj.time, obj.psth, rez, meta, cols, 'Maxdiff');
% plotStateSpace(obj.time, obj.psth, rez, cols, 'Optimization', params.dims);
lbl = {'Potent 1', 'Potent 2', 'Null 1'};
condLbl = {meta.condition{params.conditions}};
plotStateSpaceGUI(obj.time, obj.psth, rez, cols, 'Maxdiff', params.dims, lbl, condLbl);
cond = params.conditions;
plotSingleTrialsGUI(obj.time,obj.trialpsth,rez,cols,'Maxdiff', params.dims,lbl,condLbl, meta.trialNum,cond);

%% METHOD 2 ( for each cell, project onto each component, subtract from psth)

% psth = obj.psth(:,:,params.conditions);
% 
% preppsth = psth(prepix,:,:);
% preppsth = [preppsth(:,:,1) ; preppsth(:,:,2)]; % (ct,n)
% 
% [pcs,~,~,~,explained] = pca(preppsth);
% 
% % find number of pcs to explain 95% of variance
% numPCs = numComponentsToExplainVariance(explained, 90);
% rez.Qnull = pcs(:,1:numPCs); 
% 
% % subtract out each rez.Qnull latent from psth of each cell
% psth_minusprepmodes = psth;
% for i = 1:size(psth,2) % for each cell
%     for j = 1:size(psth,3) % for each condition
%         latent = sum(psth(:,i,j) * rez.Qnull(i,:), 2);
%         psth_minusprepmodes(:,i,j) = psth_minusprepmodes(:,i,j) - latent;
%     end
% end
% 
% % find potent modes
% move_minusprepmodes = psth_minusprepmodes(moveix,:,:);
% move_minusprepmodes = [move_minusprepmodes(:,:,1) ; move_minusprepmodes(:,:,2)]; % (ct,n)
% 
% [pcs,~,~,~,explained] = pca(move_minusprepmodes);
% % find number of pcs to explain 95% of variance
% numPCs = numComponentsToExplainVariance(explained, 90);
% 
% rez.Qpotent = pcs(:,1:numPCs);
% 
% % plots
% cols = {[0,0,1],[1,0,0]};
% plotLatents(obj.time, obj.psth, rez, meta, cols, 'Optimization');
% 
% dimsToPlot.potent = [1, 2];
% dimsToPlot.null = [1];
% plotStateSpace(obj.time, obj.psth, rez, 'Optimization', dimsToPlot);

end % subspaceid_maxdiff

%% 




