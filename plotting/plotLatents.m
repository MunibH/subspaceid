function plotLatents(time, psth, rez, meta, cols, methodName)

if strcmpi(methodName,'optimization')
    ylims = [-1,1];
elseif strcmpi(methodName,'maxdiff')
    ylims = [-2.2,2.2];
end

% Prep Dims
plotSubspaceLatents(time, psth, rez.Qnull, rez.Qnull_ve, meta, cols, ylims, methodName, 'Null')

% Move Dims
plotSubspaceLatents(time, psth, rez.Qpotent, rez.Qpotent_ve, meta, cols, ylims, methodName, 'Potent')


end % plotLatents

%%

function plotSubspaceLatents(time, psth, mode, varexp, meta, cols, ylims, methodName, subspace)

if size(mode,2) > 6
    ncols = 3;
else
    ncols = 2;
end

f = figure;
set(gcf, 'Position', [290   127   932   671])
for i = 1:size(mode,2) % dimension
    h(i) = subplot(ceil(size(mode,2)/2),ncols,i);
    title(['Dim ' num2str(i) '   |   %VE: ' num2str(varexp(i))]);
    xlim([meta.tmin, meta.tmax]);
    ylim(ylims);
    hold on
    for j = 1:size(psth,3) % condition
        proj = psth(:,:,j) * mode(:,i);
        plot(time, mySmooth(proj,100), 'Color', cols{j}, ...
            'LineWidth', 2.5);
    end
    hold off
end
totalVE = sum(varexp);
titleString = sprintf('%s Dimensions - %s   |   %%VE: %.2f', subspace,methodName,totalVE);
sgtitle(titleString, 'FontSize', 15)
h = axes(f, 'Visible', 'off');
h.XLabel.Visible='on';
h.YLabel.Visible='on';
xlabel('Time from go cue (s)', 'FontSize', 15)
ylabel(h,'Activity (a.u.)', 'FontSize', 15);

end % plotSubspaceLatents










