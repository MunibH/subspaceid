function load_and_process()
% load Erin's data from 2afc experiments
% output file contains alm and irn data, pooled across number of specified
% sessions/animals
% smooth, normalize, mean-center the data
% get covariance matrices of preparatory and movement epochs
clear,clc,close all

%% load data
data_pth = '/Users/Munib/Documents/Economo-Lab/subspace-id/data/2afc/ALM/';
data_file = 'data_structure_EEL6_2020-03-01.mat';
save_pth = '/Users/Munib/Documents/Economo-Lab/subspace-id/data/kaufman/';
whos_data = 'erin';

load(fullfile(data_pth,data_file));

%% parse neural data
% only use 'good' trials and clusters
trials =  get_consistent_trials(obj); % {left,right}
clusters = get_clusters_to_use(obj); 

dt = 1/400; % to match video data
t_0 = 0;
t_f = mode(obj.bp.ev.goCue)+2;
full_time = t_0:dt:t_f;
    
% get psths
probe = 1;
psth_alm = get_raw_psth(obj,probe,full_time,trials,clusters{probe});
probe = 2;
psth_irn = get_raw_psth(obj,probe,full_time,trials,clusters{probe});

% concatenate cell type psths
psth = [psth_alm ; psth_irn];
alm_idx = 1:numel(psth_alm);
irn_idx = (numel(psth_alm)+1):numel(psth);

% align time to go cue
full_time = full_time - mode(obj.bp.ev.goCue);

% % get psths for all cells, separated by trial type, and smooth
% tag.psth -> psth for all cells, all trial types
% trial types -> 1 ='hit&L&~lickearly&~stim', 2='hit&R&~lickearly&~stim'
trialTypes = {1,2};

smooth_window = 100;
% trial-avg psths of all cells in format {left psth, right psth}
% psth = (time,cells)
psth = get_psth(psth,full_time,trialTypes,smooth_window);
% remove any cells (columns) that have nans or infs (shouldn't be any)
psth = remove_inf_nan_cells(psth);

% mean firing rate for each cell
mean_fr = mean([mean(psth{1});mean(psth{2})])';

% % soft-normalize 
% firing rate of a neuron, x in (t,1), is transformed as:
% x_norm = x / (lambda + max(x) - min(x))
lambda = 0.1 * 0.005; % spike rate times delta t
psth = cellfun(@(x) x./(lambda+max(x)-min(x)),psth,'UniformOutput',false);

% regular z-scoring
% psth = cellfun(@(x) normalize(x),psth,'UniformOutput',false);

% % mean-center
% compute mean activity of each neuron across conditions at each time point
% subtract mean from each condition's response
% psth = cellfun(@(x) x-mean(x), psth,'UniformOutput',false); % mean-center within conditions
psth = mean_center(psth); % mean center across conditions
    
% % define epochs
prepInterval = 1.0; % defined as number of seconds before go cue to go cue
moveInterval = 1.0; % defined as go cue to number of seconds after go cue
[prepIdx,moveIdx] = get_epochs(full_time,prepInterval,moveInterval,dt);

% combine data
% Nprep = [left_psth(prepIdx,:);right_psth(prepIdx,:)] (ct,n)
% Nmove = [left_psth(moveIdx,:);right_psth(moveIdx,:)]
% N = {Nprep,Nmove}
[N,time] = split_data(psth,prepIdx,moveIdx,dt);

% % plot data after processing

t1 = 1:floor(length(time)/2); % first half of time vector
t2 = (floor(length(time)/2)+1):length(time); % second half

N_left = [N{1}(t1,:);N{2}(t1,:)];
N_right = [N{1}(t2,:);N{2}(t2,:)];

y_lim(1) = min([min(min(N_left)) min(min(N_right))]);
y_lim(2) = max([max(max(N_left)) max(max(N_right))]);

% figure
% subplot(1,2,1);
% plot(time,N_left); xlim([-1,1]); ylim(y_lim); title('Left')
% xlabel('time relative to go cue')
% subplot(1,2,2);
% plot(time,N_right); xlim([-1,1]); ylim(y_lim); title('Right')
% xlabel('time relative to go cue')
% sgtitle('psth - centered')

% % store vars in struct

data.psth = psth; % {left,right}
data.psth_epoch = N; % {prep,move}
data.time = full_time;
data.mean_fr = mean_fr; 
data.alm_idx = alm_idx;
data.irn_idx = irn_idx;

clearvars -except obj data data_file save_pth whos_data dt time full_time trials

%% parse movement data
% see video_features.m for video data information

% % setup views, features, coordinates, trials
views = [1,2];
viewName = {'side','bottom'};

feats = {[1,2,3],[1,2,3,4,5,6,8]};
featName = {{'tongue','jaw','nose'},...
            {'top_tongue','topleft_tongue','bottom_tongue',...
             'leftbottom_tongue','top_paw','bottom_paw','jaw'}};
coords = {[1,2],[1,2]}; % {[x,z],[x,y]}
coordName = {{'x','z'},{'x','y'}};

% which trials to use for analysis
trials = trials; % {left,right}

% % get trial-averaged feature trajectories
% trajs are normalized and mean-centered 
% data.feat_name_view = {{left_coord1,right_coord1},{left_coord2,right_coord2}}
for v = 1:numel(views) % for every view
    for f = 1:numel(feats{v}) % for every feature
        for c = 1:numel(coords{v})
            traj.([featName{v}{f} '_' viewName{v} '_' coordName{v}{c}]) = ...
                feat_traj(obj,trials,full_time,views(v),feats{v}(f),coords{v}(c));
        end
    end
end

% % create a muscle cell array, M. We'll call it feat_mat
% feat_mat = {(time,muscles_left),(time,muscles_right)}
feat_mat = cell(1,2);
feat_mat{1} = [];
feat_mat{2} = [];
feat_mat_names = fieldnames(traj);
for v = 1:numel(views) % for every view
    for f = 1:numel(feats{v}) % for every feature
        for c = 1:numel(coords{v})
            for t = 1:numel(trials)
                feat_mat{t} = [feat_mat{t} traj.([featName{v}{f} '_' viewName{v} '_' coordName{v}{c}]){t}];
            end
        end
    end
end

% % define epochs
prepInterval = 1.0; % defined as number of seconds before go cue to go cue
moveInterval = 1.0; % defined as go cue to number of seconds after go cue
[prepIdx,moveIdx] = get_epochs(full_time,prepInterval,moveInterval,dt);

% % combine data
% Mprep = [left_feat(prepIdx,:);right_feat(prepIdx,:)] (ct,n)
% Mmove = [left_feat(moveIdx,:);right_feat(moveIdx,:)]
% N = {Nprep,Nmove}
[M,time] = split_data(feat_mat,prepIdx,moveIdx,dt);

% % plot data after processing
t1 = 1:floor(length(time)/2); % first half of time vector
t2 = (floor(length(time)/2)+1):length(time); % second half

M_left = [M{1}(t1,:);M{2}(t1,:)];
M_right = [M{1}(t2,:);M{2}(t2,:)];

y_lim(1) = min([min(min(M_left)) min(min(M_right))]);
y_lim(2) = max([max(max(M_left)) max(max(M_right))]);

% figure
% subplot(1,2,1);
% plot(time,M_left); xlim([-1,1]); ylim(y_lim); title('Left')
% xlabel('time relative to go cue')
% subplot(1,2,2);
% plot(time,M_right); xlim([-1,1]); ylim(y_lim); title('Right')
% xlabel('time relative to go cue')
% sgtitle('feature matrix - centered')

% % store vars in struct

data.feat_mat = feat_mat; % {left,right}
data.feat_names = feat_mat_names;
data.feat_mat_epoch = M; % {prep,move}
data.prep_idx = prepIdx;
data.move_idx = moveIdx;
data.trial_types = {'left','right'};
data.epoch_types = {'prep','move'};

clearvars -except data data_file save_pth whos_data

%% Reduce dimensionality of N(psth) and M(feat_mat)
% N = 6 dims, M = 3 dims    ( for now )
d_N = 6;
d_M = 3;

% % perform pca on N = [left_N,right_N] (t x cn)
% perform pca on M = [left_M,right_M] (t x cn)
[~,N_red] = pca([data.psth{1}(:,data.irn_idx),data.psth{2}(:,data.irn_idx)]);
N_red = N_red(:,1:d_N);
[~,M_red] = pca([data.feat_mat{1},data.feat_mat{2}]);
M_red = M_red(:,1:d_M);

% % plot N_red/M_red
% figure
% subplot(2,1,1); plot(data.time,N_red); title('N')
% subplot(2,1,2); plot(data.time,M_red); title('M')

data.N = N_red; % use prep_idx and move_idx to get Nprep,move...
data.M = M_red;

%% save

data.parent = data_file;
save_fn = [data.parent(16:end-4) '_processed_irn.mat'];
save(fullfile(save_pth,whos_data,save_fn),'data')

end % loadAndProcess

%% Helper Functions

function psth = get_psth(psth,time,tt,smooth_window)
% get psths for every cell, separated by trial type
    left_psth = zeros(numel(time),numel(psth));
    right_psth = left_psth;
    for cell = 1:numel(psth) % for every cell
        left_psth(:,cell) = my_smooth(mean(psth{cell}{tt{1}},2),smooth_window);
        right_psth(:,cell) = my_smooth(mean(psth{cell}{tt{2}},2),smooth_window);
    end
    psth = {left_psth,right_psth};
end % get_psth


function psth_norm = my_normalize(psth,tag)
% firing rate of a neuron, x in (t,1), is normalized as:
% x_norm = (x - mean(x(baselineIdx)) / std(x(baselineIdx))
    % get baseline idx from tag
    baselineIdx = tag.epoch.ix{1};
    psth_norm = cell(1,numel(psth));
    for tt = 1:numel(psth) % for each trial type
        psth_norm{tt} = zeros(size(psth{tt}));
        for cellIdx = 1:size(psth{tt},2) % for every cell
            mu = mean(psth{tt}(baselineIdx,cellIdx));
            stdev = std(psth{tt}(baselineIdx,cellIdx));
            psth_norm{tt}(:,cellIdx) = (psth{tt}(:,cellIdx) - mu) / stdev;
        end
    end
end % my_normalize

function [prepIdx,moveIdx] = get_epochs_full(tag,prepInterval,moveInterval,dt)
% get the time index of the prep and move epochs
    prepInterval = prepInterval / dt; 
    moveInterval = moveInterval / dt; 
    % find index of go cue
    [~,t0_idx] = min(abs(tag.time-0));
    % define epochs relative to go cue index
    prepIdx = floor(t0_idx-prepInterval:t0_idx);
    moveIdx = floor(t0_idx:t0_idx+moveInterval);
end % get_epochs_full

function [prepIdx,moveIdx] = get_epochs(time,prepInterval,moveInterval,dt)
% get the time index of the prep and move epochs
    prepInterval = prepInterval / dt; 
    moveInterval = moveInterval / dt; 
    % find index of go cue
    [~,t0_idx] = min(abs(time-0));
    % define epochs relative to go cue index
    prepIdx = floor(t0_idx-30-prepInterval:t0_idx-30);
    moveIdx = floor(t0_idx+1:t0_idx+1+moveInterval);
end % get_epochs

function [N,time] = split_data(psth,prepIdx,moveIdx,dt)
% combine data such that 
% Nprep = [left_psth(prepIdx,:);right_psth(prepIdx,:)] (ct,n)
% Nmove = [left_psth(moveIdx,:);right_psth(moveIdx,:)]
% N = {Nprep,Nmove}
    N = cell(1,2);
    N{1} = [psth{1}(prepIdx,:);psth{2}(prepIdx,:)]; % prep
    N{2} = [psth{1}(moveIdx,:);psth{2}(moveIdx,:)]; % move
% return a new time vector
    time = -(length(prepIdx)*dt):dt:0;
    time = [time dt:dt:((length(moveIdx)*dt)-dt)];
end % split_data

function psth_scrubbed = remove_inf_nan_cells(psth)
    psth_scrubbed = cell(1,numel(psth));
    for tt = 1:numel(psth)
        [~,c_nan] = find(isnan(psth{tt}));
        [~,c_inf] = find(isinf(psth{tt}));
        to_delete = unique([c_nan,c_inf]);
        psth{tt}(:,to_delete) = [];
        psth_scrubbed{tt} = psth{tt};
    end
end % remove_inf_nan_cells

function plot_cov(data,cutout)
    num_to_cut = ceil( numel(data) * cutout / 2);
    sorted_data = sort(data(:));
    cmin = sorted_data( num_to_cut );
    cmax = sorted_data( end - num_to_cut + 1);
    imagesc(data, [cmin, cmax]);
    colorbar;
end

function centered = mean_center(psth)
% compute mean activity of each neuron across conditions at each time point
% subtract mean from each condition's response
    centered = psth;
    for cluIdx = 1:size(psth{1},2)
        for t = 1:size(psth{1},1)
            % mean across conditions at time t
            mean_cond_t = mean( [psth{1}(t,cluIdx) psth{2}(t,cluIdx)] );
            centered{1}(t,cluIdx) = centered{1}(t,cluIdx) - mean_cond_t;
            centered{2}(t,cluIdx) = centered{2}(t,cluIdx) - mean_cond_t;
        end
    end
end

function psth = get_raw_psth(obj,probe,time,trials,clusters)
    numClusters = length(clusters);
    psth = cell(numClusters, 1); % {left, right}
    for cluIdx = 1:numClusters
        psth{cluIdx} = cell(1,numel(trials)); % {left,right})
        data = obj.clu{probe}(clusters(cluIdx));
        for tt = 1:numel(trials)
            psth{cluIdx}{tt} = zeros(length(time),numel(trials{tt}));
            for trialIdx = 1:numel(trials{tt})
                trialtm = data.trialtm(ismember(data.trial, trials{tt}(trialIdx)));
                psth{cluIdx}{tt}(:,trialIdx) = histc(trialtm,time);
            end
        end
    end
end % get_raw_psth

function trialsToUse = get_consistent_trials(obj)
% narrow down the trials to use to trials where delay and go cue occur at
% most frequent time and where there were successful licks. 
    % get trials where delay epoch starts at mode of obj.bp.ev.delay
    % and go cue starts at mode of obj.bp.ev.goCue
    delay = obj.bp.ev.delay; % (nTrials,1) double 
    goCue = obj.bp.ev.goCue; % (nTrials,1) double 
    modeDelay = mode(delay);
    modeGoCue = mode(goCue);
    trialsMask = (delay==modeDelay) & (goCue==modeGoCue);    
    % get successful lick trials
    right_hit_trials = obj.bp.hit & obj.bp.R;
    left_hit_trials = obj.bp.hit & obj.bp.L;
    trialsMask_right = trialsMask & right_hit_trials;
    trialsMask_left = trialsMask & left_hit_trials;    
    % get trials to use
    trials_right = find(trialsMask_right);
    trials_left = find(trialsMask_left);        
    % {right,left}
    trialsToUse = {trials_left, trials_right};
end % get_consistent_trials

function clustersToUse = get_clusters_to_use(obj)
% get clusters of quality not equal to 'Poor' to use for analysis
    clustersToUse = {zeros(1,500),zeros(1,500)};
    for probe = 1:numel(obj.clu)
        numClusters = numel(obj.clu{probe});
        ct = 1;
        for cluster = 1:numClusters
            if ~strcmp(obj.clu{1}(cluster).quality,'Poor')
                clustersToUse{probe}(ct) = cluster;
                ct = ct + 1;
            end
        end
        % remove unused space
        clustersToUse{probe}(clustersToUse{probe}==0) = [];    
    end
end % get_clusters_to_use

function traj = feat_traj(obj,trials,time,view,feat,coord)
% get trial-by-trial feature trajectories
    traj = cell(1,2); % {left_coord,right_coord}
    for cond = 1:numel(trials)
        traj{cond} = zeros(length(time),length(trials{cond}));
        for trialIdx = 1:length(trials{cond})
            % get trial data
            traj{cond}(:,trialIdx) = ...
                obj.traj{view}(trials{cond}(trialIdx)).ts(1:length(time),coord,feat);
            % fill missing data
            traj{cond}(:,trialIdx) = fillmissing(traj{cond}(:,trialIdx),'nearest');
            % normalize b/w 0 and 1
            traj{cond}(:,trialIdx) = (traj{cond}(:,trialIdx) - mean(traj{cond}(:,trialIdx)))...
                / (max(traj{cond}(:,trialIdx)) - min(traj{cond}(:,trialIdx)));
            % smooth
            traj{cond}(:,trialIdx) = my_smooth(traj{cond}(:,trialIdx),10);
        end
    end
% trial-average
    for cond = 1:numel(trials)
        traj{cond} = mean(traj{cond},2);
    end
    
end % feat_traj



