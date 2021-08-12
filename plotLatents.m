function plotLatents(time, psth, rez, meta, cols, methodName)

% Prep Dims
ylims = [-3,3];
plotSubspaceLatents(time, psth, rez.Qnull, rez.Qnull_ve, meta, cols, ylims, methodName, 'Null')

% Move Dims
ylims = [-3,3];
plotSubspaceLatents(time, psth, rez.Qpotent, rez.Qpotent_ve, meta, cols, ylims, methodName, 'Potent')


end % plotLatents

%%

function plotSubspaceLatents(time, psth, mode, varexp, meta, cols, ylims, methodName, subspace)

f = figure;
set(gcf, 'Position', [290   127   932   671])
for i = 1:size(mode,2)
    h(i) = subplot(ceil(size(mode,2)/2),2,i);
    title(['Dim ' num2str(i) '   |  %VE: ' num2str(varexp(i))]);
    xlim([meta.tmin, meta.tmax]);
    ylim(ylims);
    hold on
    for j = 1:size(psth,3)
        proj = psth(:,:,j) * mode(:,i);
        plot(time, proj, 'Color', cols{j}, ...
            'LineWidth', 2.5);
    end
    hold off
end
sgtitle([subspace ' Dimensions - ' methodName])
h = axes(f, 'Visible', 'off');
h.XLabel.Visible='on';
h.YLabel.Visible='on';
xlabel('Time from go cue (s)', 'FontSize', 20)
ylabel(h,'Activity (a.u.)', 'FontSize', 20);

end % plotSubspaceLatents










