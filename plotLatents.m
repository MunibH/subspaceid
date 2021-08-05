function plotLatents(time, psth, rez, meta, cols, methodName)

% Prep Dims
figure;
for i = 1:size(rez.Qnull,2)
    subplot(size(rez.Qnull,2),1,i)
    xlim([meta.tmin, meta.tmax]);
    hold on
    for j = 1:size(psth,3)
        proj = psth(:,:,j) * rez.Qnull(:,i);
        plot(time, proj, 'Color', cols{j});
    end
    hold off
end
sgtitle(['Null Dimensions - ' methodName])
xlabel('Time from go cue (s)')

% Move Dims
figure;
for i = 1:size(rez.Qpotent,2)
    subplot(size(rez.Qpotent,2),1,i)
    xlim([meta.tmin, meta.tmax]);
    hold on
    for j = 1:size(psth,3)
        proj = psth(:,:,j) * rez.Qpotent(:,i);
        plot(time, proj, 'Color', cols{j});
    end
    hold off
end
sgtitle(['Potent Dimensions - ' methodName])
xlabel('Time from go cue (s)')

end % plotLatents



