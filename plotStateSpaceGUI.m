function plotStateSpaceGUI(time, psth, rez, cols, methodName, dims)

%% SETUP GUI
bcol = [1 1 1];

% Main Figure
h.fig(1) = figure(1000);
set(h.fig(1), 'Units', 'Pixels', 'Position', [828 379 1032 751], 'Color', bcol);

h.popupmenu(i) = uicontrol('Style', 'popupmenu', 'Units', 'pixels', 'Position', ...
        [375 775-35*i 100 25], 'String', h.feat.str,'Callback', ...
        {@updateVideo,gcf}, 'Value', i, 'BackgroundColor', [1 1 1]);

end %plotStateSpaceGUI