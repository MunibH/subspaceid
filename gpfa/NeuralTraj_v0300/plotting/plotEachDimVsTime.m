function plotEachDimVsTime(seq, xspec, time, varargin)
%
% plotEachDimVsTime(seq, xspec, binWidth, ...)
%
% Plot each state dimension versus time in a separate panel.
%
% INPUTS:
%
% seq       - data structure containing extracted trajectories
% xspec     - field name of trajectories in 'seq' to be plotted 
%             (e.g., 'xorth' or 'xsm')
% binWidth  - spike bin width used when fitting model
%
% OPTIONAL ARGUMENTS:
%
% nPlotMax  - maximum number of trials to plot (default: 20)
% redTrials - vector of trialIds whose trajectories are plotted in red
%             (default: [])
% nCols     - number of subplot columns (default: 4)
%
% @ 2009 Byron Yu -- byronyu@stanford.edu

  nPlotMax  = 20;
  redTrials = [];
  nCols     = 4;
  assignopts(who, varargin);

  f = figure;
  pos = get(gcf, 'position');
  set(f, 'position', [pos(1) pos(2)./2 2*pos(3) pos(4)]);

  Xall = [seq.(xspec)];
  xMax = ceil(10 * max(abs(Xall(:)))) / 10; % round max value to next highest 1e-1
  
  Tmax    = max([seq.T]);  

  ytk     = [-xMax 0 xMax];

  nRows   = ceil(size(Xall, 1) / nCols);
  
  for n = 1:min(length(seq), nPlotMax)
    dat = seq(n).(xspec);
    T   = seq(n).T;
        
    for k = 1:size(dat,1)
      subplot(nRows, nCols, k);
      hold on;
      
      if ismember(seq(n).trialId, redTrials)
        col = [1 0 0]; % red
        lw  = 1;
      else
        col = 0.2 * [1 1 1]; % gray
        lw = 1;
      end      
      plot(time, dat(k,:), 'linewidth', lw, 'color', col);
    end
  end

  for k = 1:size(dat,1)
    h = subplot(nRows, nCols, k);
    axis([min(time) max(time) 1.1*min(ytk) 1.1*max(ytk)]);

    if isequal(xspec, 'xorth')
      str = sprintf('$$\\tilde{\\mathbf x}_{%d,:}$$',k);
    else
      str = sprintf('$${\\mathbf x}_{%d,:}$$',k);
    end
    title(str, 'interpreter', 'latex', 'fontsize', 16);
        
    set(h, 'ytick', ytk, 'yticklabel', ytk);
    xlabel('Time (s) from Go Cue');
  end
