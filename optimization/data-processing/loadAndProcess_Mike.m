function loadAndProcess_Mike()
% load mike economo's data from Nature 2018 paper : Distinct descending motor cortex pathways and their roles in movement
% contains 1207 cells (tagged + untagged), but the two populations are
% unidentified
% smooth, normalize, mean-center the data
% get covariance matrices of preparatory and movement epochs
clear,clc,close all

%% load data
data_pth = '/Users/Munib/Documents/Economo-Lab/subspace-id/data/elsayed/';
whos_data = 'mike';
data_file = 'alldat.mat'; % all data (includes tagged cells but unidentified)
% data_file = 'alldat_tagged.mat'; % all data including manually merged tagged cells
load(fullfile(data_pth,whos_data,data_file));

%% get psths for all cells, separated by trial type, and smooth
% tag.psth -> psth for all cells, all trial types
% trial types -> 1 ='hit&L&~lickearly&~stim', 2='hit&R&~lickearly&~stim'
trialTypes = {1,2};

smooth_window = 25;
% trial-avg psths of all cells in format {left psth, right psth}
% psth = (time,cells)
psth = get_psth(tag,trialTypes,smooth_window);
% remove any cells (columns) that have nans or infs (shouldn't be any)
psth = remove_inf_nan_cells(psth);

psth_tavg = psth;

% save normalized firing rates for each cell for subspace contribution analysis
mean_fr = mean([mean(psth{1});mean(psth{2})])';

% %% filter to make less squiggly
% 
% Fs = 30000;
% N = size(psth{1},1);
% Nfir = 70;
% Fst = 75;
% firf = designfilt('lowpassfir','FilterOrder',Nfir, ...
%     'CutoffFrequency',Fst,'SampleRate',Fs);
% psthf = cellfun(@(x) filter(firf,x),psth,'UniformOutput',false);
% 
% % correct time lag
% delay = mean(grpdelay(firf));
% 
% tt = tag.time(1:end-delay);
% psthn = cellfun(@(x) x(1:end-delay,:),psth,'UniformOutput',false);
% psthff = psthf;
% psthff{1}(1:delay,:) = [];
% psthff{2}(1:delay,:) = [];
% psth = psthff;
% tag.time = tt;

%% soft-normalize 
% firing rate of a neuron, x in (t,1), is transformed as:
% x_norm = x / (lambda + max(x) - min(x))
lambda = 10 * 0.005; % spike rate times delta t
psth = cellfun(@(x) x./(lambda+max(x)-min(x)),psth,'UniformOutput',false);
% regular z-scoring
% psth = cellfun(@(x) normalize(x),psth,'UniformOutput',false);
psth_sn = psth;
% save off XL and XR to project on dims
psth_int = cellfun(@(x) x-mean(x), psth,'UniformOutput',false);

%% mean-center
% compute mean activity of each neuron across conditions at each time point
% subtract mean from each condition's response
psth = mean_center(psth); % mean-center across conditions

% psth = cellfun(@(x) x-mean(x), psth,'UniformOutput',false); % mean-center within conditions


%% define epochs
prepInterval = 1.0; % defined as number of seconds before go cue to go cue
moveInterval = 1.0; % defined as go cue to number of seconds after go cue
dt = mode(diff(tag.time)); % time step
[prepIdx,moveIdx] = get_epochs(tag,prepInterval,moveInterval,dt);

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
Cprep = cov(N{1});
Cmove = cov(N{2});

cutout = 0.01;

figure
subplot(2,2,1)
plot_cov(Cprep,cutout); title('Prep Covariance')
subplot(2,2,2)
plot_cov(Cmove,cutout); title('Move Covariance')
subplot(2,2,3)
imagesc(imresize(Cprep,1/6,'bilinear')); title('Prep Downsampled'); colorbar;
subplot(2,2,4)
imagesc(imresize(Cmove,1/6,'bilinear')); title('Move Downsampled'); colorbar;

%% PT Lower and PT Upper cells
if strcmp(data_file,'alldat_tagged.mat')
    % get psth and covariance matrices of ptl and ptu cells
    psth_ptlow = {psth{1}(:,tag.ptlow_ix),psth{2}(:,tag.ptlow_ix)};
    psth_ptup = {psth{1}(:,tag.ptup_ix),psth{2}(:,tag.ptup_ix)};
    [N_ptlow,~] = split_data(psth_ptlow,prepIdx,moveIdx,dt);
    [N_ptup,~] = split_data(psth_ptup,prepIdx,moveIdx,dt);
    Cprep_ptlow = cov(N_ptlow{1});
    Cmove_ptlow = cov(N_ptlow{2});
    Cprep_ptup = cov(N_ptup{1});
    Cmove_ptup = cov(N_ptup{2});
end

%% save

data.full_psth = psth;
data.psth_int = psth_int;
data.psth_sn = psth_sn;
data.psth_tavg = psth_tavg;
data.full_time = tag.time;
data.time = time;
data.psth = N;
data.mean_fr = mean_fr; 
data.Cprep = Cprep;
data.Cmove = Cmove;
data.pop_idx = 1:1207; % first 1207 cells (not tagged)

if strcmp(data_file,'alldat_tagged.mat')
    data.ptlow_idx = tag.ptlow_ix;
    data.ptup_idx = tag.ptup_ix;
    data.psth_ptlow = psth_ptlow;
    data.psth_ptup = psth_ptup;
    data.Cprep_ptlow = Cprep_ptlow;
    data.Cmove_ptlow = Cmove_ptlow;
    data.Cprep_ptup = Cprep_ptup;
    data.Cmove_ptup = Cmove_ptup;
    save_fn = 'alldat_tagged_for_grant.mat'; % 1337 cells (1207 + ptlower + ptupper)
else
    save_fn = 'alldat_processed.mat'; % 1207 cells
end

save(fullfile(data_pth,whos_data,save_fn),'data')

end % loadAndProcess

%% Helper Functions

function psth = get_psth(tag,tt,smooth_window)
% get psths for every cell, separated by trial type
    left_psth = zeros(numel(tag.time),numel(tag.psth));
    right_psth = left_psth;
    for cell = 1:numel(tag.psth) % for every cell
        left_psth(:,cell) = mySmooth(mean(tag.psth{cell}{tt{1}},2),smooth_window);
        right_psth(:,cell) = mySmooth(mean(tag.psth{cell}{tt{2}},2),smooth_window);
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

function [prepIdx,moveIdx] = get_epochs(tag,prepInterval,moveInterval,dt)
% get the time index of the prep and move epochs
    prepInterval = prepInterval / dt; 
    moveInterval = moveInterval / dt; 
    % find index of go cue
    [~,t0_idx] = min(abs(tag.time-0));
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
end % plot_cov

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
end % mean_center















