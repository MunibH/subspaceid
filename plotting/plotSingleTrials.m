function plotSingleTrials(fig, time, psth, rez, cols, methodName, dims, lbl, condlbl, trials, condnum)

h = guidata(fig);

% can only plot in 3D
assert((numel(dims.potent) + numel(dims.null)) == 3, 'Plot 3 dimensions!')

% compute latents
nullproj = nan(numel(time), size(psth,3), size(rez.Qnull,2));
potentproj = nan(numel(time), size(psth,3), size(rez.Qpotent,2));
for i = 1:size(rez.trialpsth,3) % for each trial
    % (time,trials,dims)
    nullproj(:,i,1:rez.dPrep) = (rez.trialpsth(:,:,i) * rez.Qnull);
    potentproj(:,i,1:rez.dMove) = (rez.trialpsth(:,:,i) * rez.Qpotent);
end

% get specific potent and null dims to plot
if ~isempty(dims.potent)
    dat.potent = potentproj(:,:,dims.potent);
end
if ~isempty(dims.null)
    dat.null = nullproj(:,:,dims.null);
end

fn = fieldnames(dat);
ct = 1;
for i = 1:numel(fn)
    for j = 1:size(dat.(fn{i}),3)
        dat.toPlot(:,:,ct) = dat.(fn{i})(:,:,j);
        ct = ct + 1;
    end
end


% plot
nTrialsToPlot = 20; % per condition
smth = str2double(h.smooth.String);
for i = 1:numel(condnum)
    trix = trials{condnum(i)};
    trix = datasample(trix,nTrialsToPlot,'Replace',false);
    
    p(1:nTrialsToPlot,i) = ...
           plot3(mySmooth(dat.toPlot(:,trix,1),smth), ...
                 mySmooth(dat.toPlot(:,trix,2),smth), ...
                 mySmooth(dat.toPlot(:,trix,3),smth), ... 
                 'LineWidth', 0.5, ...
                 'DisplayName', condlbl{i}, ...
                 'Color', cols{i});
         
    hold on
    
end
hold off
grid on
xlabel(lbl{1}, 'FontSize', 30)
ylabel(lbl{2}, 'FontSize', 30)
zlabel(lbl{3}, 'FontSize', 30)
title(methodName, 'FontSize', 30);
legend(p(1,:), ...
       'Location', 'bestoutside','Units','normalized', ...
       'Position', [0.703716285422975,0.804930894683553,0.262596892871598,0.092543273092745])

% create gradient color trajectories
stretch = str2double(h.colstretch.String);
alph = 0.2;
righttrajcol = [uint8(summer(size(psth,1)*stretch)*255) , ...
                uint8(alph*ones(size(psth,1)*stretch,1))].';
lefttrajcol = [uint8(autumn(size(psth,1)*stretch)*255) , ...
                uint8(alph*ones(size(psth,1)*stretch,1))].';
drawnow
for i = 1:numel(condnum)
    if i == 1
        trajcol = righttrajcol;
    elseif i == 2
        trajcol = lefttrajcol;
    end
    for j = 1:nTrialsToPlot
        set(p(j,i).Edge, 'ColorBinding','interpolated', 'ColorData',trajcol(:,1:size(psth,1)))
    end
end




end