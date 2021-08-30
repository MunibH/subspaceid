function rez = subspaceid_maxdiff_singletrials(meta, obj, params)

%% PREP and MOVE EPOCHS
% prepix and moveix corresponds to time idx for each epoch
[rez.prepix, rez.moveix] = getEpochs(obj.time, meta.prepEpoch, meta.moveEpoch);

%% get null and potent data
nulldat = cell(obj.bp.Ntrials,1);
potentdat = cell(obj.bp.Ntrials,1);
allix = 1:numel(obj.time);
nanmat = false(size(obj.trialpsth));
for i = 1:obj.bp.Ntrials
    movix = meta.movix{i};
    dat = mySmooth(obj.trialpsth(:,:,i), 51);
    potentdat{i} = dat(movix,:); % all clusters, 1 trial
    nonmovix = find(~ismember(allix,movix));
    nulldat{i} = dat(nonmovix,:); % all clusters, 1 trial
    nanmat(movix, :, i) = true;
end

null = zeros(size(nanmat, 1), size(nanmat, 2), numel(meta.trialNum));
for i = 1:numel(meta.trialNum)
    dat =  obj.trialpsth(:, :, meta.trialNum{i});
    dat(nanmat(:,:,meta.trialNum{i})) = NaN;
    null(:, :, i) = nanmean(dat, 3);
end
null = null(rez.prepix, :, :);
null = reshape(permute(null, [1 3 2]), size(null, 1)*numel(meta.trialNum), 72);
null = null-mean(null, 1);

potent = zeros(size(nanmat, 1), size(nanmat, 2), numel(meta.trialNum));
for i = 1:numel(meta.trialNum)
    dat =  obj.trialpsth(:, :, meta.trialNum{i});
    potent(:, :, i) = nanmean(dat, 3);
end
potent = potent(rez.moveix,:,:);
potent = reshape(permute(potent, [1 3 2]), size(potent, 1)*numel(meta.trialNum), 72);
potent = potent - mean(potent, 1);


%% pca on null data, then potent data, then orthogonalize potent to null
nullmodes = pca(null,'NumComponents',params.xDim);

nulllatents = nan(numel(obj.time),params.xDim,obj.bp.Ntrials);
for i = 1:obj.bp.Ntrials
    nulllatents(:,:,i) = (obj.trialpsth(:,:,i)) * nullmodes;
end
alph = 0.4;
ncols = 3;
lw = 1;
figure; hold on;
for i = 1:params.xDim
    subplot(ceil(size(nullmodes,2)/3),ncols,i); hold on
%     plot(obj.time, mySmooth(squeeze(nulllatents(:,i,meta.trialNum{1})), 51),'Color',[0,0,1,alph]);
%     plot(obj.time, mySmooth(squeeze(nulllatents(:,i,meta.trialNum{2})), 51),'Color',[1,0,0,alph]);
%     plot(obj.time, mySmooth(squeeze(nulllatents(:,i,meta.trialNum{3})), 51),'Color',[0,1,1,alph]);
%     plot(obj.time, mySmooth(squeeze(nulllatents(:,i,meta.trialNum{4})), 51),'Color',[1,0,1,alph]);
    
    plot(obj.time, mean(squeeze(nulllatents(:,i,meta.trialNum{1})), 2),'b', 'LineWidth', lw)
    plot(obj.time, mean(squeeze(nulllatents(:,i,meta.trialNum{2})), 2),'r', 'LineWidth', lw)
    plot(obj.time, mean(squeeze(nulllatents(:,i,meta.trialNum{3})), 2),'c', 'LineWidth', lw)
    plot(obj.time, mean(squeeze(nulllatents(:,i,meta.trialNum{4})), 2),'m', 'LineWidth', lw)

end
sgtitle('null modes')

% pca on potent data and then orthogonalize
potentmodes = pca(potent,'NumComponents',params.xDim);
% orthogonalize to null modes
orthModes = gschmidt([nullmodes potentmodes]); % where mode is (clusters,dims)
potentmodes = orthModes(:,(params.xDim+1):end);
 
potentlatents = nan(numel(obj.time),params.xDim,obj.bp.Ntrials);
for i = 1:obj.bp.Ntrials
    potentlatents(:,:,i) = obj.trialpsth(:,:,i) * potentmodes;
end
figure; hold on;
for i = 1:params.xDim
    subplot(ceil(size(potentmodes,2)/3),ncols,i); hold on
%     plot(obj.time, mySmooth(squeeze(potentlatents(:,i,meta.trialNum{1})), 25),'Color',[0,0,1,alph]);
%     plot(obj.time, mySmooth(squeeze(potentlatents(:,i,meta.trialNum{2})), 25),'Color',[1,0,0,alph]);
    plot(obj.time, mean(squeeze(potentlatents(:,i,meta.trialNum{1})), 2),'b', 'LineWidth', lw)
    plot(obj.time, mean(squeeze(potentlatents(:,i,meta.trialNum{2})), 2),'r', 'LineWidth', lw)
    plot(obj.time, mean(squeeze(potentlatents(:,i,meta.trialNum{3})), 2),'c', 'LineWidth', lw)
    plot(obj.time, mean(squeeze(potentlatents(:,i,meta.trialNum{4})), 2),'m', 'LineWidth', lw)
end
sgtitle('potent modes orthog to null')


%% proj single trials onto null modes, subtract the proj from single trials
% % then do pca on the residual to get potent mode
% removeNull = nan([size(obj.trialpsth) params.xDim]); % (time,clu,trials,xdims)
% for cluix = 1:size(nullmodes,1)
%     for dim = 1:params.xDim
%         removeNull(:,cluix,:,dim) = squeeze(obj.trialpsth(:,cluix,:)) * nullmodes(cluix,dim);
%     end
% end
% removeNull = sum(removeNull,4);
% 
% potent = obj.trialpsth - removeNull;
% % concatenate in time
% potent = permute(potent,[1 3 2]);
% potent = reshape(potent,size(potent,1)*size(potent,2),size(potent,3));
% 
% % pca on potent data
% potentmodes = pca(potent,'NumComponents',params.xDim);
%  
% potentlatents = nan(numel(obj.time),params.xDim,obj.bp.Ntrials);
% for i = 1:obj.bp.Ntrials
%     potentlatents(:,:,i) = obj.trialpsth(:,:,i) * potentmodes;
% end
% figure; hold on;
% for i = 1:params.xDim
%     subplot(ceil(size(potentmodes,2)/3),ncols,i); hold on
%     plot(obj.time, mySmooth(squeeze(potentlatents(:,i,meta.trialNum{1})), 25),'Color',[0,0,1,alph]);
%     plot(obj.time, mySmooth(squeeze(potentlatents(:,i,meta.trialNum{2})), 25),'Color',[1,0,0,alph]);
%     plot(obj.time, mean(squeeze(potentlatents(:,i,meta.trialNum{2})), 2),'m', 'LineWidth', 3)
%     plot(obj.time, mean(squeeze(potentlatents(:,i,meta.trialNum{1})), 2),'c', 'LineWidth', 3)
% end
% sgtitle('potent modes subtracting null')




%% pca on potent data, then pca on null data, then orthogonalize null to potent
potentmodes = pca(potent,'NumComponents',params.xDim);

 potentlatents = nan(numel(obj.time),params.xDim,obj.bp.Ntrials);
for i = 1:obj.bp.Ntrials
    potentlatents(:,:,i) = obj.trialpsth(:,:,i) * potentmodes;
end
figure; hold on;
for i = 1:params.xDim
    subplot(ceil(size(potentmodes,2)/3),ncols,i); hold on
%     plot(obj.time, mySmooth(squeeze(potentlatents(:,i,meta.trialNum{1})), 25),'Color',[0,0,1,alph]);
%     plot(obj.time, mySmooth(squeeze(potentlatents(:,i,meta.trialNum{2})), 25),'Color',[1,0,0,alph]);
    plot(obj.time, mean(squeeze(potentlatents(:,i,meta.trialNum{1})), 2),'b', 'LineWidth', lw)
    plot(obj.time, mean(squeeze(potentlatents(:,i,meta.trialNum{2})), 2),'r', 'LineWidth', lw)
    plot(obj.time, mean(squeeze(potentlatents(:,i,meta.trialNum{3})), 2),'c', 'LineWidth', lw)
    plot(obj.time, mean(squeeze(potentlatents(:,i,meta.trialNum{4})), 2),'m', 'LineWidth', lw)
end
sgtitle('potent modes')

% pca on null data
nullmodes = pca(null,'NumComponents',params.xDim);
% orthogonalize to null modes
orthModes = gschmidt([potentmodes nullmodes]); % where mode is (clusters,dims)
nullmodes = orthModes(:,(params.xDim+1):end);

nulllatents = nan(numel(obj.time),params.xDim,obj.bp.Ntrials);
for i = 1:obj.bp.Ntrials
    nulllatents(:,:,i) = (obj.trialpsth(:,:,i)) * nullmodes;
end

figure; hold on;
for i = 1:params.xDim
    subplot(ceil(size(nullmodes,2)/3),ncols,i); hold on
%     plot(obj.time, mySmooth(squeeze(nulllatents(:,i,meta.trialNum{1})), 51),'Color',[0,0,1,alph]);
%     plot(obj.time, mySmooth(squeeze(nulllatents(:,i,meta.trialNum{2})), 51),'Color',[1,0,0,alph]);
%     plot(obj.time, mySmooth(squeeze(nulllatents(:,i,meta.trialNum{3})), 51),'Color',[0,1,1,alph]);
%     plot(obj.time, mySmooth(squeeze(nulllatents(:,i,meta.trialNum{4})), 51),'Color',[1,0,1,alph]);
    
    plot(obj.time, mean(squeeze(nulllatents(:,i,meta.trialNum{1})), 2),'b', 'LineWidth', lw)
    plot(obj.time, mean(squeeze(nulllatents(:,i,meta.trialNum{2})), 2),'r', 'LineWidth', lw)
    plot(obj.time, mean(squeeze(nulllatents(:,i,meta.trialNum{3})), 2),'c', 'LineWidth', lw)
    plot(obj.time, mean(squeeze(nulllatents(:,i,meta.trialNum{4})), 2),'m', 'LineWidth', lw)

end
sgtitle('null modes orthog to potent')

end % subspaceid_maxdiff_singletrials

%% 




