clear,clc,close all
data_pth = 'Y:\EEL\Experiments\EEL17\Analysis';
date = '2021-03-10';
load(fullfile(data_pth,date,'Trajectories.mat'));

med = 7;
view = 2;
    
ntrials = 254;
len = 5000;  
ts1 = NaN.*zeros(len, 3, ntrials);

trials = ones(1,ntrials);

trialnum = find(trials);
for i = 1:ntrials
    temp = traj{view}(trialnum(i)).ts(:, :, 1)+traj{view}(trialnum(i)).ts(:, :, 3);
    dat = temp/2;
    ts1(1:size(dat, 1), :, i) = dat;
end


ts = medfilt1(ts1, med, [], 1);

cols = {(1/255)*[0,0,255],(1/255)*[7,138,195],(1/255)*[0,255,0],...
        (1/255)*[255,0,0],(1/255)*[155,40,40]};

% figure;

% get every lick trajectory
for trix = 1:size(ts,3) % for every trial
    trialdat = ts(:,1,trix);
    y = isnan(trialdat);
    
    licknum = 1;
    curlickix = 1;
    for i = 1:length(y)
        if y(i) == 1
            continue
        else
            idx(licknum,curlickix) = i;
            curlickix = curlickix + 1;
            if i ~=length(y)
                if y(i+1) == 1
                    licknum = licknum + 1;
                end
            end
        end
    end

    % convert idx to a cell array
    % only keep first nlicks (first lick only for block position)
    nlicks = size(idx,1);
    for i = 1:nlicks
        dat = idx(i,:);
        cellidx(i) = {dat(dat~=0)};
    end

    % plot first lick for each trial 
    for i = 1:numel(cellidx)
        numel(cellidx)
        numel(cellidx{1})
        'hi'
        plot(ts(cellidx{i},1,trix),-ts(cellidx{i},2,trix));
%         plot(ts(cellidx{i},1,trix),-ts(cellidx{i},2,trix),'-','Color',cols{i});
        hold on
    end
    'hi'
end
hold off



