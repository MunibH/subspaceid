function plotStateSpace(h, time, psth, rez, cols, methodName, dims, lbl)

% can only plot in 3D
assert((numel(dims.potent) + numel(dims.null)) == 3, 'Plot 3 dimensions!')

% compute latents
null_latent = nan(numel(time), size(rez.Qnull,2), size(psth,3));
potent_latent = nan(numel(time), size(rez.Qpotent,2), size(psth,3));
for i = 1:size(psth,3)
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
[~,gocueidx] = min(abs(time-0));
for i = 1:size(dat.toPlot,3)
    p(i) = plot3(dat.toPlot(:,1,i), dat.toPlot(:,2,i),...
             dat.toPlot(:,3,i), 'LineWidth', 2); hold on
    plot3(dat.toPlot(gocueidx,1,i), dat.toPlot(gocueidx,2,i),...
             dat.toPlot(gocueidx,3,i), '.', 'Color', cols{i}, 'MarkerSize',50);
end
hold off
grid on
xlabel(lbl{1}, 'FontSize', 30)
ylabel(lbl{2}, 'FontSize', 30)
zlabel(lbl{3}, 'FontSize', 30)
title(methodName, 'FontSize', 30);

% create gradient color trajectories
stretch = 1.5;
righttrajcol = [uint8(winter(size(psth,1)*stretch)*255) , ...
                uint8(ones(size(psth,1)*stretch,1))].';
lefttrajcol = [uint8(hot(size(psth,1)*stretch)*255) , ...
                uint8(ones(size(psth,1)*stretch,1))].';
drawnow
set(p(1).Edge, 'ColorBinding','interpolated', 'ColorData',righttrajcol(:,1:size(psth,1)))
set(p(2).Edge, 'ColorBinding','interpolated', 'ColorData',lefttrajcol(:,1:1:size(psth,1)))


end % plotStateSpace



