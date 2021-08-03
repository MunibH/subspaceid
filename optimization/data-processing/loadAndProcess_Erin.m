function loadAndProcess_Erin()
% load Erin's data from 2afc experiments
% output file contains alm and irn data, pooled across number of specified
% sessions/animals
% smooth, normalize, mean-center the data
% get covariance matrices of preparatory and movement epochs
clear,clc,close all

%% load data
data_pth = '/Users/Munib/Documents/Economo-Lab/subspace-id/data/2afc/ALM/';
data_file = {'data_structure_EEL6_2020-03-01.mat','data_structure_EEL7_2020-03-02.mat'};
save_pth = '/Users/Munib/Documents/Economo-Lab/subspace-id/elsayed_v2/data/';
whos_data = 'erin';

combined_alm = cell(1,numel(data_file));
combined_irn = cell(1,numel(data_file));
for i = 1:numel(data_file)
    load(fullfile(data_pth,data_file{i}));

    %% get stuff
    trials =  get_consistent_trials(obj); % {left,right}
    clusters = get_clusters_to_use(obj); 

    dt = 0.005; % same as mike's data
    t_0 = 0;
    t_f = mode(obj.bp.ev.goCue)+2;
    full_time = t_0:dt:t_f;
    
    % get psths of all cells
    probe = 1;
    combined_alm{i} = get_raw_psth(obj,probe,full_time,trials,clusters{probe});
    probe = 2;
    combined_irn{i} = get_raw_psth(obj,probe,full_time,trials,clusters{probe});

end
% concatenate combined_psths
psth_alm = {};
psth_irn = {};
for i = 1:numel(data_file)
    psth_alm = [psth_alm;combined_alm{i}];
    psth_irn = [psth_irn;combined_irn{i}];
end
psth = [psth_alm;psth_irn];
alm_idx = 1:numel(psth_alm);
irn_idx = ((numel(psth_alm)+1)):numel(psth);

% align time to go cue
full_time = full_time - mode(obj.bp.ev.goCue);

%% get psths for all cells, separated by trial type, and smooth
% tag.psth -> psth for all cells, all trial types
% trial types -> 1 ='hit&L&~lickearly&~stim', 2='hit&R&~lickearly&~stim'
trialTypes = {1,2};

smooth_window = 100;
% trial-avg psths of all cells in format {left psth, right psth}
% psth = (time,cells)
psth = get_psth(psth,full_time,trialTypes,smooth_window);
% remove any cells (columns) that have nans or infs (shouldn't be any)
psth = remove_inf_nan_cells(psth);

% save normalized firing rates for each cell for subspace contribution analysis
mean_fr = mean([mean(psth{1});mean(psth{2})])';

%% soft-normalize 
% firing rate of a neuron, x in (t,1), is transformed as:
% x_norm = x / (lambda + max(x) - min(x))
lambda = 10 * 0.005; % spike rate times delta t
psth = cellfun(@(x) x./(lambda+max(x)-min(x)),psth,'UniformOutput',false);

% regular z-scoring
% psth = cellfun(@(x) normalize(x),psth,'UniformOutput',false);

%% mean-center
% compute mean activity of each neuron across conditions at each time point
% subtract mean from each condition's response
% psth = cellfun(@(x) x-mean(x), psth,'UniformOutput',false); % mean-center within conditions
psth = mean_center(psth); % mean center across conditions
    
%% define epochs
prepInterval = 1.0; % defined as number of seconds before go cue to go cue
moveInterval = 1.0; % defined as go cue to number of seconds after go cue
[prepIdx,moveIdx] = get_epochs(full_time,prepInterval,moveInterval,dt);

% combine data
% Nprep = [left_psth(prepIdx,:);right_psth(prepIdx,:)] (ct,n)
% Nmove = [left_psth(moveIdx,:);right_psth(moveIdx,:)]
% N = {Nprep,Nmove}
[N,time] = split_data(psth,prepIdx,moveIdx,dt);

t1 = 1:floor(length(time)/2); % first half of time vector
t2 = (floor(length(time)/2)+1):length(time); % second half

N_left = [N{1}(t1,:);N{2}(t1,:)];
N_right = [N{1}(t2,:);N{2}(t2,:)];

%% plot data after processing

y_lim(1) = min([min(min(N_left)) min(min(N_right))]);
y_lim(2) = max([max(max(N_left)) max(max(N_right))]);

figure
subplot(1,2,1);
plot(time,N_left); xlim([-1,1]); ylim(y_lim); title('Left')
xlabel('time relative to go cue')
subplot(1,2,2);
plot(time,N_right); xlim([-1,1]); ylim(y_lim); title('Right')
xlabel('time relative to go cue')
sgtitle('psth - centered')

%% get covariance matrices
Cprep = cov(N{1}); % all cells (alm + irn)
Cmove = cov(N{2});
Cprep_alm = cov(N{1}(:,alm_idx)); % alm cells
Cmove_alm = cov(N{2}(:,alm_idx));
Cprep_irn = cov(N{1}(:,irn_idx)); % irn cells
Cmove_irn = cov(N{2}(:,irn_idx));

% cutout = 0.01;
% 
% figure
% subplot(2,2,1)
% plot_cov(Cprep,cutout); title('Prep Covariance')
% subplot(2,2,2)
% plot_cov(Cmove,cutout); title('Move Covariance')
% subplot(2,2,3)
% imagesc(imresize(Cprep,1/6,'bilinear')); title('Prep Downsampled'); colorbar;
% subplot(2,2,4)
% imagesc(imresize(Cmove,1/6,'bilinear')); title('Move Downsampled'); colorbar;


%% save

data.full_psth = psth;
data.full_time = full_time;
data.time = time;
data.psth = N; % {left,right}
data.mean_fr = mean_fr; 
data.psth_alm = {psth{1}(:,alm_idx),psth{2}(:,alm_idx)};
data.psth_irn = {psth{1}(:,irn_idx),psth{2}(:,irn_idx)};
data.Cprep = Cprep;
data.Cmove = Cmove;
data.Cprep_alm = Cprep_alm;
data.Cmove_alm = Cmove_alm;
data.Cprep_irn = Cprep_irn;
data.Cmove_irn = Cmove_irn;
data.alm_idx = alm_idx;
data.irn_idx = irn_idx;

save_fn = 'EEL62020-03-01_EEL72020-03-02_processed_softnorm.mat'; % description necessary

save(fullfile(save_pth,whos_data,save_fn),'data')

end % loadAndProcess

%% Helper Functions

function psth = get_psth(psth,time,tt,smooth_window)
% get psths for every cell, separated by trial type
    left_psth = zeros(numel(time),numel(psth));
    right_psth = left_psth;
    for cell = 1:numel(psth) % for every cell
        left_psth(:,cell) = mySmooth(mean(psth{cell}{tt{1}},2),smooth_window);
        right_psth(:,cell) = mySmooth(mean(psth{cell}{tt{2}},2),smooth_window);
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



