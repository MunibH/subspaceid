function subspaceid_maxdiff(meta, obj, params)

%% PREPROCESS
[obj, meta] = preprocess_maxdiff(meta, obj);

%% PREP and MOVE EPOCHS
% prepix and moveix corresponds to time idx for each epoch
[prepix, moveix] = getEpochs(obj.time, meta.prepEpoch, meta.moveEpoch);

%% METHOD 1 ( N*(I-WW') )

%% FIND NULL MODES
psth = obj.psth(:,:,params.conditions);

preppsth = psth(prepix,:,:);
preppsth = [preppsth(:,:,1) ; preppsth(:,:,2)]; % (ct,n)

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
moveproj = proj(moveix,:,:);
moveproj = [moveproj(:,:,1) ; moveproj(:,:,2)]; % (ct,n)

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
    xlim([meta.tmin, meta.tmax]);
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
    xlim([meta.tmin, meta.tmax]);
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
    p(i) = plot3(potent_latent(:,1,i), potent_latent(:,2,i),...
             null_latent(:,1,i), 'LineWidth', 2);
end
hold off
grid on
xlabel('Potent 1')
ylabel('Potent 2')
zlabel('Null 1')
% modified jet-colormap
stretch = 1.5;
righttrajcol = [uint8(winter(size(psth,1)*stretch)*255) uint8(ones(size(psth,1)*stretch,1))].';
lefttrajcol = [uint8(hot(size(psth,1)*stretch)*255) uint8(ones(size(psth,1)*stretch,1))].';
drawnow
set(p(1).Edge, 'ColorBinding','interpolated', 'ColorData',righttrajcol(:,1:size(psth,1)))
set(p(2).Edge, 'ColorBinding','interpolated', 'ColorData',lefttrajcol(:,1:1:size(psth,1)))

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
% 
% % Prep Dims
% figure;
% for i = 1:size(rez.Qnull,2)
%     subplot(size(rez.Qnull,2),1,i)
%     hold on
%     for j = 1:size(psth,3)
%         proj = psth(:,:,j) * rez.Qnull(:,i);
%         plot(obj.time, proj, 'Color', cols{j});
%     end
%     hold off
% end
% sgtitle('Null Dimensions')
% xlabel('Time from go cue (s)')
% 
% % Move Dims
% figure;
% for i = 1:size(rez.Qpotent,2)
%     subplot(size(rez.Qpotent,2),1,i)
%     hold on
%     for j = 1:size(psth,3)
%         proj = psth(:,:,j) * rez.Qpotent(:,i);
%         plot(obj.time, proj, 'Color', cols{j});
%     end
%     hold off
% end
% sgtitle('Potent Dimensions')
% xlabel('Time from go cue (s)')
% 
% % State Space (Potents dims 1 and 2 and null dim 1)
% null_latent = nan(numel(obj.time), size(rez.Qnull,2), size(psth,3));
% potent_latent = nan(numel(obj.time), size(rez.Qpotent,2), size(psth,3));
% for i = 1:size(psth,3)
%     null_latent(:,:,i) = psth(:,:,i) * rez.Qnull; 
%     potent_latent(:,:,i) = psth(:,:,i) * rez.Qpotent; 
% end
% 
% figure; hold on
% for i = 1:size(psth,3)
%     plot3(potent_latent(:,1,i), potent_latent(:,2,i),...
%            null_latent(:,1,i), 'Color', cols{i})
% end
% hold off
% xlabel('Potent 1')
% ylabel('Potent 2')
% zlabel('Null 1')

end % subspaceid_maxdiff

%% 




