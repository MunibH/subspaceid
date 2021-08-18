function mov = getMoveIdx(obj)
% returns cell array mov where each entry is an (Nframesx1) logical vector
% labeling where in a trial the animal is moving

% algorithm:
% First, set timepoint to true if tongue is visible
% Then, set timepoint to true if jaw speed exceeds a threshold
% Last, set timepoint true if variance in jaw position exceeds a threshold

%% params

dt = 1/400;
view = 1; % side cam
tongfeat = 1;
jawfeat = 2; 
sm = 7; % smoothing window


vid = obj.traj{view};
mov = cell(obj.bp.Ntrials,1);
for j = 1:obj.bp.Ntrials
    %% get tongue trajectory
    
    p = vid(j).ts(:, 3, tongfeat);
    
    Nframes = size(vid(j).ts, 1);
    dat = zeros(Nframes, 2);
    dat(:,1) = (1:Nframes)./(1/dt); % time
    dat(:,2) = mySmooth(vid(j).ts(:,1,tongfeat), sm); % x
    dat(:,3) = mySmooth(vid(j).ts(:,2,tongfeat), sm); % z
    dat(p<0.9,[2,3]) = NaN;
    
    %% find where tongue is visible 
    
    tongmov = false(Nframes,1);
    tongmov(~isnan(dat(:,2))) = true;
    tongmov(~isnan(dat(:,3))) = true;
    
    %% get jaw trajectory

    p = vid(j).ts(:, 3, jawfeat);
    
    Nframes = size(vid(j).ts, 1);
    dat = zeros(Nframes, 2);
    dat(:,1) = (1:Nframes)./(1/dt); % time
    dat(:,2) = mySmooth(vid(j).ts(:,1,jawfeat), sm); % x
    dat(:,3) = mySmooth(vid(j).ts(:,2,jawfeat), sm); % y
    dat(p<0.9,[2,3]) = NaN;
    
    %% fill nans
    
    for i = 2:3
        dat(:,i) = fillmissing(dat(:,i),'makima');
    end
    
    %% mean-center and normalize jaw x and z pos
    
    for i = 2:3
        dat(:,i) = (dat(:,i) - nanmin(dat(:,i))) / (nanmax(dat(:,i)) - nanmin(dat(:,i)));
        dat(:,i) = dat(:,i) - nanmean(dat(:,i));
        dat(1:5,i) = dat(6,i);
    end    
    
    
    %% get jaw speed
    
    vx = gradient(dat(:,2), dt);
    vy = gradient(dat(:,3), dt);
    dat(:,4) = sqrt(vx.^2 + vy.^2);
%     dat(:,4) = gradient(dat(:,3), dt);

    dat(1:50,4) = dat(51,4); 
    dat(:,4) = mySmooth(abs(dat(:,4)),200);
    
    dat(:,4) = (dat(:,4) - nanmin(dat(:,4))) / (nanmax(dat(:,4)) - nanmin(dat(:,4)));
    dat(:,4) = dat(:,4) - nanmean(dat(:,4));
    dat(1:50,4) = dat(51,4); 

    %% find where jaw speed exceeds a threshold
    
    per1 = 0.025;
    if min(dat(:,4)) < 0
        jawspeedthresh = per1*(max(dat(:,4)) + min(dat(:,4)));
    else
        jawspeedthresh = per1*(max(dat(:,4)) - min(dat(:,4)));
    end
    jawspeedmov = dat(:,4) > jawspeedthresh;
        
    
    %% find where variance in jaw z pos exceeds a threshold
    
    val = zeros(Nframes, 1);
    binsize = 0.5; % in seconds
    for ii = (1/dt * binsize):Nframes % for each half second bin
        val(ii) = nanvar(dat(ii-(1/dt*binsize-1):ii,3)); % variance in jaw z pos in bin
    end
    
    per2 = 0.05;
    jawzthresh = per2*max(val);
    jawzmov = val > jawzthresh;
    
    %% get mov
    
    mov{j} = tongmov | jawspeedmov | jawzmov;

    
    %% plot
%     close all
%     figure('units','normalized','Position',[0.2 0.1 0.7 0.8]); 
%     subplot(3,1,1)
%     plot(dat(:,1),dat(:,4)); hold on;
%     plot(dat(jawspeedmov,1),dat(jawspeedmov,4),'.'); hold off
%     yline(jawspeedthresh,'r');
%     ylabel('jaw speed')
%     
%     subplot(3,1,2)
%     title('variance in jaw z pos')
%     plot(dat(:,1),val); hold on
%     plot(dat(jawzmov,1),val(jawzmov),'.'); hold off
%     yline(jawzthresh,'r');
%     ylabel('variance in jaw z pos')
%     
%         subplot(3,1,3)
%     plot(dat(:,1),dat(:,3)); hold on
%     plot(dat(mov,1),dat(mov,3),'.'); hold off;
%     ylabel('jaw z position')
%     xlabel('time')
%     sgtitle(['Trial ' num2str(j)])
%     pause

end

end % getMoveidx










