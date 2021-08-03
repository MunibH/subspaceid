function subspaceid_maxdiff(meta, obj, params)

%% PREPROCESS
[obj, meta] = preprocess_maxdiff(meta, obj);

%% PREP and MOVE EPOCHS
% prepix and moveix corresponds to time idx for each epoch
[prepix, moveix] = getEpochs(obj.time, meta.prepEpoch, meta.moveEpoch);

%% FIND NULL MODES
psth = obj.psth(:,:,params.conditions);

preppsth = mean(psth(prepix,:,:),3);

[pcs,~,~,~,explained] = pca(preppsth);

% find number of pcs to explain 95% of variance
numPCs = numComponentsToExplainVariance(explained, 90);
rez.Qnull = pcs(:,1:numPCs); 

%% PROJECT OUT PROMINENT NULL MODES
modesToKeep = eye(size(pcs,1)) - (rez.Qnull*rez.Qnull');

proj = nan(size(psth));
for i = 1:size(psth,3)
    proj(:,:,i) = psth(:,:,i) * modesToKeep;
end

%% FIND POTENT MODES AS PCs OF MOVE EPOCH
moveproj = mean(proj(moveix,:,:),3);
[pcs,~,~,~,explained] = pca(moveproj);
% find number of pcs to explain 95% of variance
numPCs = numComponentsToExplainVariance(explained, 90);

rez.Qpotent = pcs(:,1:numPCs);

%% PLOTS

cols = {[0,0,1],[1,0,0]};

% Prep Dims
figure;
for i = 1:size(rez.Qnull,2)
    subplot(size(rez.Qnull,2),1,i)
    hold on
    for j = 1:size(psth,3)
        proj = psth(:,:,j) * rez.Qnull(:,i);
        plot(obj.time, proj, 'Color', cols{j});
    end
    hold off
end
sgtitle('Null Dimensions')
xlabel('Time from go cue (s)')

% Move Dims
figure;
for i = 1:size(rez.Qpotent,2)
    subplot(size(rez.Qpotent,2),1,i)
    hold on
    for j = 1:size(psth,3)
        proj = psth(:,:,j) * rez.Qpotent(:,i);
        plot(obj.time, proj, 'Color', cols{j});
    end
    hold off
end
sgtitle('Potent Dimensions')
xlabel('Time from go cue (s)')

% State Space (Potents dims 1 and 2 and null dim 1)
null_latent = nan(numel(obj.time), size(rez.Qnull,2), size(psth,3));
potent_latent = nan(numel(obj.time), size(rez.Qpotent,2), size(psth,3));
for i = 1:size(psth,3)
    null_latent(:,:,i) = psth(:,:,i) * rez.Qnull; 
    potent_latent(:,:,i) = psth(:,:,i) * rez.Qpotent; 
end

figure; hold on
for i = 1:size(psth,3)
    plot3(potent_latent(:,1,i), potent_latent(:,2,i),...
           null_latent(:,1,i), 'Color', cols{i})
end
hold off
xlabel('Potent 1')
ylabel('Potent 2')
zlabel('Null 1')

'hi'

end % subspaceid_maxdiff

%% 




