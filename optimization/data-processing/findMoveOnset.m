function findMoveOnset()
%findMoveOnset Uses behavior data to determine time of movement onset
%   For analysis, movement onset will be estimated offline based on jaw speed. 
%   Movement onset (for each trial) was defined as the time that jaw speed first surpassed 15% 
%   of the peak speed of first lick for that trial. 
%   Similarly, movement end was defined as the time when jaw speed first dropped 
%   below 15% of the peak jaw speed.

clear,clc,close all

%% load data

% if ispc
%     drive_pth = '';
% else
%     drive_pth = '/Volumes/eng_research_economo'; %/MAH/Experiments/MAH3/Analysis/2021-03-24'
% end
% who = 'MAH';
% animal = '3';
% session = '2021-03-24';
% 
% data_pth = fullfile(drive_pth,who,'Experiments',[who animal],'Analysis',session);
% data_file = dir(fullfile(data_pth,'data_structure*.mat'));
% 
% load(fullfile(data_pth,data_file.name))

data_pth = '/Users/Munib/Documents/Economo-Lab/subspace-id/data/elsayed/munib';
data_file = dir(fullfile(data_pth,'data_structure*.mat'));
load(fullfile(data_pth,data_file.name))


% vars
dt = 1/400;
view = 1; % side cam
feat = 2; % jaw
coords = [1,2]; % x,z coord
smth = 20; % smoothing window

% get index in obj.traj.ts where first lick ends for each trial
flview = 2; % bottom view, tongue, x coord
flfeat = 1;
flcoord = 1;
firstLickEndIdx = getFirstLickEndIdx(obj,flview,flfeat,flcoord,dt);

% get movement trajectories until end of first lick after go cue
traj = getTraj(obj,view,feat,coords,smth,firstLickEndIdx); % {(trials,time) coord1, (trials,time) coord2}
% time = 0:dt:(size(traj{1},2)*dt);
time = linspace(0,size(traj{1},2)*dt,size(traj{1},2));

% get speed of trajectories
trajSpeed = getTrajSpeed(traj,dt);

% find move onset time for each trial
moveOn = getMoveOnTime(trajSpeed,time,obj.bp.ev.goCue); % (trials,2 (on,off))

% latency distribution
hist(moveOn - obj.bp.ev.goCue)

'hi'

end % findMoveOnset


%% Helper Functions

function idx = getFirstLickEndIdx(obj,view,feat,coord,dt)
    idx = nan(obj.bp.Ntrials,1);
    for trial = 1:obj.bp.Ntrials
        % get trial data
        trialData = obj.traj{view}(trial).ts(:,coord,feat);
        % find time/index of go cue
        gocueTime = obj.bp.ev.goCue(trial);
        time = linspace(0,length(trialData)*dt,length(trialData));
        [~,gocueIdx] = min(abs(time - gocueTime));
        % find where tongue trajectory is nan
        y = isnan(trialData);
        % find where tongue is visible
        yis0 = find(y==0);
        % find all indices where tongue goes from visible to not visible
        endLickIdx = find(diff(yis0)~=1);
        
        % if endLickIdx is empty, that means there was only one lick
        if isempty(endLickIdx)
            % check that idx is greater than go cue idx
            if yis0(end) > gocueIdx
                idx(trial) = yis0(end);
                continue
            else
                continue
            end
        end
        
        % if multiple licks, get the first lick corresponding to after go cue
        idxs = yis0(endLickIdx);
        idxs(end+1) = yis0(end); % need to hanlde case where only one lick after go cue 
        idxs = idxs(idxs > gocueIdx);
        
        % there may not be any licks after the go cue
        try
            idx(trial) = idxs(1);
        catch
            continue
        end

    end
end % getFirstLickEndIdx

function traj = getTraj(obj,view,feat,coords,smth,lickIdx)
    % find longest trial
    maxTrialLen = max(lickIdx);

    % get data until first lick
    traj = cell(1,numel(coords));
    for coord = 1:numel(coords)
        traj{coord} = nan(obj.bp.Ntrials,maxTrialLen); % trials,time
        for trial = 1:obj.bp.Ntrials
            len = lickIdx(trial);
            % handle trials where mouse didn't lick (lickIdx(trial) == nan)
            try
                traj{coord}(trial,1:len) = obj.traj{view}(trial).ts(1:len,coord,feat); 
            catch
                continue
            end
        end
    end
    % normalize
    traj{1} = traj{1}';
    traj{2} = traj{2}';
    traj{1} = (traj{1} - min(traj{1})) ./ (max(traj{1}) - min(traj{1}));
    traj{2} = (traj{2} - min(traj{2})) ./ (max(traj{2}) - min(traj{2}));
    traj{1} = traj{1}';
    traj{2} = traj{2}';
    % smooth
    traj{1} = mySmooth(traj{1}',smth)';
    traj{2} = mySmooth(traj{2}',smth)';
end % getTraj


function trajSpeed = getTrajSpeed(traj,dt)
    trajSpeed = zeros(size(traj{1},1),size(traj{1},2)); %(trials,time)
    for i = 1:size(traj{1},1)
        x_vel = my_diff(traj{1}(i,:),dt);
        z_vel = my_diff(traj{2}(i,:),dt);
        trajSpeed(i,:) = sqrt(x_vel.^2 + z_vel.^2);
%         trajSpeed(i,:) = sqrt(z_vel.^2);
    end

end % trajSpeed


function moveOn = getMoveOnTime(speed,time,goCue)
    moveOn = nan(size(speed,1),1); % (trials,1)
    for trial = 1:size(speed,1)
        
        % omit trials where video cut off before trial ended (these were
        % really long trials where the mouse wouldn't stop licking)
        if goCue(trial) > time(end)
            continue
        end
        
        % for each trial, get 0.2 sec before and after after go cue
        [~,goCueIdx] = min(abs(time - (goCue(trial))));
        trialTraj = speed(trial,goCueIdx:end);
        
        % get 15% of peak speed
        thresh = 0.15*max(trialTraj);
        % find when speed after go cue exceeds thresh
        idx = find(trialTraj > thresh);
        % skip trials with no licks after go cue
        try
            moveOnIdx = idx(1) + goCueIdx - 1;
        catch
            continue
        end
        moveOn(trial) = time(moveOnIdx); % time of movement onset
        
%         plot(time,speed(trial,:),'LineWidth',2); hold on
%         xlim([goCue(trial)-0.2,goCue(trial)+0.3]);
%         xline(goCue(trial),'k','LineWidth',1); hold on
%         plot(moveOn(trial),speed(trial,moveOnIdx),'r.', ...
%             'MarkerSize',20); hold on
%         yline(thresh,'r'); hold off
%         title(['Trial ' num2str(trial)]);
%         pause
    end
    
end % getMoveOnTime




