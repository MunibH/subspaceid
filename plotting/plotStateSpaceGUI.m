function plotStateSpaceGUI(time, psth, rez, cols, methodName, dims, lbl, cond)

%% SETUP GUI
bcol = [1 1 1];

% Main Figure
h.fig(1) = figure(1000);  hold on

set(h.fig(1), 'Units', 'Pixels', 'Position', [828 379 1032 751], 'Color', bcol);

% axis menus
menuString = createMenuString(rez);
textString = {'X: ', 'Y: ', 'Z: '};
for i = 1:3
    uicontrol('Parent',h.fig(1), 'Style','text','String', textString{i}, 'Position',[820 600-75*i 30 25], 'Fontsize', 15, ...
        'BackgroundColor', bcol)
    h.axismenu(i) = uicontrol('Parent',h.fig(1), 'Style', 'popupmenu', 'Units', 'pixels', 'Position', ...
        [850 600-75*i 100 25], 'String', menuString, 'Value', i, 'BackgroundColor', bcol);
end

% set defaults
for i = 1:3
    a = h.axismenu(i);
    strToFind = {lbl{i}};

    fun = @(s)~cellfun('isempty',strfind(a.String,s));
    out = cellfun(fun,strToFind,'UniformOutput',false);
    idx = find(all(horzcat(out{:}),2));
    if numel(idx) > 1
        idx = idx(1);
    end
    set(h.axismenu(i),'Value',idx)
end

% smooth
uicontrol('Parent',h.fig(1), 'Style','text','String', 'Smooth: ', 'Position',[820 250 95 25], 'Fontsize', 15, ...
        'BackgroundColor', bcol)
h.smooth = uicontrol('Parent',h.fig(1), 'Style', 'edit', 'Position', [920 250 50 25], ...
    'String', '0', 'BackgroundColor', bcol, ...
    'Callback', {@updateDims,gcf,time,psth,rez,cols,methodName,cond});
h.smooth.String = '100';

% color stretching
uicontrol('Parent',h.fig(1), 'Style','text','String', 'ColStretch: ', 'Position',[820 200 95 25], 'Fontsize', 15, ...
        'BackgroundColor', bcol)
h.colstretch = uicontrol('Parent',h.fig(1), 'Style', 'edit', 'Position', [920 200 50 25], ...
    'String', '1', 'BackgroundColor', bcol, ...
    'Callback', {@updateDims,gcf,time,psth,rez,cols,methodName,cond});

% update plot button
h.update = uicontrol('Parent',h.fig(1), 'Style', 'pushbutton', 'Position', [815 315 180 40], ...
    'String', 'Plot Dims','FontWeight','Bold', 'BackgroundColor', [0.3, 0.7, 0.5], ...
    'Callback', {@updateDims,gcf,time,psth,rez,cols,methodName,cond});

% axes
set(gca,'Visible','off')
h.ax = axes;
set(h.ax, 'Position', [0.0835 0.0974 0.6551 0.8307]); 

guidata(h.fig(1),h)

% initial plot
plotStateSpace(h.fig(1), time, psth, rez, cols, methodName, dims, lbl, cond);


set(h.ax,'Color',[0 0 0 0.2])



end %plotStateSpaceGUI

%%

function menuString = createMenuString(rez)

menuString = cell(rez.dPrep + rez.dMove, 1);
for i = 1:rez.dPrep
    menuString(i) = {['Null ' num2str(i)]};
end
for j = 1:rez.dMove
    menuString(j+rez.dPrep) = {['Potent ' num2str(j)]};
end

end % createMenuString


function updateDims(~,~,fig,time,psth,rez,cols,methodName,cond)

h = guidata(fig);

dims.potent = [];
dims.null = [];

dimString = cell(numel(h.axismenu),1);
for i = 1:numel(h.axismenu)
    dimString(i) = h.axismenu(i).String(h.axismenu(i).Value);
    dimStringSplit = strsplit(dimString{i});
    num = str2double(dimStringSplit{2});
    if contains(dimString{i}, 'Potent')
        dims.potent = [dims.potent num];
    elseif contains(dimString{i}, 'Null')
        dims.null = [dims.null num];
    end
    
end

plotStateSpace(fig, time, psth, rez, cols, methodName, dims, dimString, cond);

set(h.ax,'Color',[0 0 0 0.2])

end % updateDims

















