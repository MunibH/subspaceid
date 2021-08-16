function plotStateSpace(fig, time, psth, rez, cols, methodName, dims, lbl, cond)

h = guidata(fig);

% can only plot in 3D
assert((numel(dims.potent) + numel(dims.null)) == 3, 'Plot 3 dimensions!')

% compute latents
null_latent = nan(numel(time), size(rez.Qnull,2), size(psth,3));
potent_latent = nan(numel(time), size(rez.Qpotent,2), size(psth,3));
for i = 1:size(psth,3)
    % (time,dims,numconds)
    null_latent(:,:,i) = psth(:,:,i) * rez.Qnull; 
    potent_latent(:,:,i) = psth(:,:,i) * rez.Qpotent; 
end

% get specific potent and null dims to plot
if ~isempty(dims.potent)
    dat.potent = potent_latent(:,dims.potent,:);
end
if ~isempty(dims.null)
    dat.null = null_latent(:,dims.null,:);
end

fn = fieldnames(dat);
ct = 1;
for i = 1:numel(fn)
    for j = 1:size(dat.(fn{i}),2)
        dat.toPlot(:,ct,:) = dat.(fn{i})(:,j,:);
        ct = ct + 1;
    end
end

% plot
smth = str2double(h.smooth.String);
[~,gocueidx] = min(abs(time-0));
for i = 1:size(dat.toPlot,3)
    p(i) = plot3(mySmooth(dat.toPlot(:,1,i),smth), ...
                 mySmooth(dat.toPlot(:,2,i),smth), ...
                 mySmooth(dat.toPlot(:,3,i),smth), ... 
                 'LineWidth', 2, ...
                 'DisplayName', cond{i});
         
    hold on
    
    plot3(dat.toPlot(gocueidx,1,i), dat.toPlot(gocueidx,2,i),...
             dat.toPlot(gocueidx,3,i), '.', 'Color', cols{i}, ...
             'DisplayName', 'Go Cue', 'MarkerSize',50);
end
hold off
grid on
xlabel(lbl{1}, 'FontSize', 30)
ylabel(lbl{2}, 'FontSize', 30)
zlabel(lbl{3}, 'FontSize', 30)
title(methodName, 'FontSize', 30);
legend('Location', 'bestoutside','Units','normalized', ...
       'Position', [0.703716285422975,0.804930894683553,0.262596892871598,0.092543273092745])

% create gradient color trajectories
stretch = str2double(h.colstretch.String);
righttrajcol = [uint8(winter(size(psth,1)*stretch)*255) , ...
                uint8(ones(size(psth,1)*stretch,1))].';
lefttrajcol = [uint8(hot(size(psth,1)*stretch)*255) , ...
                uint8(ones(size(psth,1)*stretch,1))].';
drawnow
set(p(1).Edge, 'ColorBinding','interpolated', 'ColorData',righttrajcol(:,1:size(psth,1)))
set(p(2).Edge, 'ColorBinding','interpolated', 'ColorData',lefttrajcol(:,1:1:size(psth,1)))


end % plotStateSpace



