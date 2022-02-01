function rez = subspaceid_maxdiff(meta, obj, params)

%% PREPROCESS
rez.psth = obj.psth(:,:,params.condToUse);
rez.time = obj.time;
rez.condition = params.condition(params.condToUse);
rez.trialid = params.trialid(params.condToUse);

rez = preprocess_maxdiff(rez);

%% PREP and MOVE EPOCHS
% prepix and moveix corresponds to time idx for each epoch
[rez.prepix, rez.moveix] = getEpochs(obj.time, params.prepEpoch, params.moveEpoch);

%% METHOD 1 ( N*(I-WW') )

%% FIND NULL MODES

preppsth = rez.psth(rez.prepix,:,:);
preppsth = [preppsth(:,:,1) ; preppsth(:,:,2)]; % (ct,n)

% pca
[pcs,lambda,explained] = myPCA(preppsth);

% find number of pcs to explain 95% of variance
rez.dPrep = numComponentsToExplainVariance(explained, params.varToExplain);
rez.Qnull = pcs(:,1:rez.dPrep);
rez.Qnull_ve = explained(1:rez.dPrep); % var explained of prep epoch in null subspace

%% PROJECT OUT PROMINENT NULL MODES
modesToKeep = eye(size(pcs,1)) - (rez.Qnull*rez.Qnull');

proj = nan(size(rez.psth));
for i = 1:size(rez.psth,3)
    proj(:,:,i) = rez.psth(:,:,i) * modesToKeep;
end

%% FIND POTENT MODES AS PCs OF MOVE EPOCH
% moveproj = proj(rez.moveix,:,:); % only use move epoch for potent mode ID
moveproj = proj; % use all leftover data for potent mode ID (seems more right)
moveproj = [moveproj(:,:,1) ; moveproj(:,:,2)]; % (ct,n)

% pca
[pcs,lambda,explained] = myPCA(moveproj);

% find number of pcs to explain 95% of variance of moveproj (not psth)
rez.dMove = numComponentsToExplainVariance(explained, params.varToExplain);
rez.Qpotent = pcs(:,1:rez.dMove);

% variance explained of move epoch by potent space
movepsth = rez.psth(rez.moveix,:,:);
movepsth = [movepsth(:,:,1) ; movepsth(:,:,2)]; % (ct,n)
[~,lambdafull,~] = myPCA(movepsth);
explained = (lambda / sum(lambdafull)) * 100;
rez.Qpotent_ve = explained(1:rez.dMove);

%% PLOTS

% cols = params.cols;
% plotLatents(obj.time, obj.psth, rez, meta, cols, 'Maxdiff');
% lbl = {'Potent 1', 'Potent 2', 'Null 1'};
% condLbl = {meta.condition{params.conditions}};
% plotStateSpaceGUI(obj.time, obj.psth, rez, cols, 'Maxdiff', params.dims, lbl, condLbl);
% cond = params.conditions;
% plotSingleTrialsGUI(obj.time,obj.trialpsth,rez,cols,'Maxdiff', params.dims,lbl,condLbl, meta.trialNum,cond);


clrs = {[0 0 1],[1 0 0],[0.5 0.5 1],[1 0.5 0.5],'k'};

Q = rez.Qnull;

for i = 1:size(Q,2) % dimension
    f = figure;
    title(['Dim ' num2str(i) '   |   %VE: ' num2str(rez.Qnull_ve(i))]);
    xlim([params.tmin, params.tmax]);
    hold on
    for j = 1:size(rez.psth,3) % condition
        proj = rez.psth(:,:,j) * Q(:,i);
        plot(obj.time, mySmooth(proj,100), 'Color', clrs{j}, ...
            'LineWidth', 2.5);
    end
    hold off
    ax = gca;
    ax.FontSize = 20;
end

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




