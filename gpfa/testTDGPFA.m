function testTDGPFA()
close all;
addpath(genpath('C:\Users\labadmin\Desktop\gpfa\NeuralTraj_v0300'));
parent = 'C:\Users\labadmin\Desktop\gpfa';
fn = 'data_structure_JEB7_2021-04-29.mat';
temp = load(fullfile(parent, fn));
obj = temp.obj;
clear temp

dt = 0.005;
probenum = 1;
edges = 0:dt:5;
midpts = edges(1:end-1) + mean(diff(edges))./2;
Ntpts = numel(midpts);
Nclust = numel(obj.clu{probenum});

if contains(fn,'JEB7')
    tr = find(obj.bp.hit&~obj.bp.early&(obj.bp.autowater.nums'==2)); % JEB7
else
    tr = find(obj.bp.hit&~obj.bp.early&~obj.bp.autowater); % EKH2
end
Ntrials = numel(tr);

% % Remove low-firing rate units, e.g., all those firing less than 5
%   spikes per second on average across all trials.
%   
%   The fitted observation noise (diagonal element of R) for a
%   low-firing rate unit will be small, so the neural trajectory may
%   show a deflection each time the neuron spikes.
psth = zeros(Ntpts,Nclust);
for i = 1:Nclust
    spkix = ismember(obj.clu{probenum}(i).trial, tr);
    if numel(spkix)==0 || numel(spkix)==1  % kilisort outputted clusters with no spikes??
        continue
    end
    N = histc(obj.clu{probenum}(i).trialtm(spkix), edges);
    N = N(1:end-1);
    if size(N,1) < size(N,2)
        N = N';
    end
    psth(:,i) = MySmooth(N./numel(tr)./dt, 15);

end
meanFRs = mean(psth);
lowFR = 5;
use = meanFRs > lowFR;

% badClusts = [26 57 62 92]; % JEB7
% badClusts = [85 109 112 133 174 189]; % MAH8
% use(badClusts) = false;

obj.clu{probenum} = obj.clu{probenum}(use);
Nclust = numel(obj.clu{probenum});

% format data for gpfa
for i = 1:Ntrials
    dat(i).trialId = tr(i);
    dat(i).spikes = zeros(Nclust, Ntpts);
    for j = 1:Nclust
        N = histc(obj.clu{probenum}(j).trialtm(obj.clu{probenum}(j).trial==tr(i)), edges);
        N = N(1:end-1);
        dat(i).spikes(j, :) = N; %clusters by time points
    end
end

% GPFA
runIdx = 2; 
method = 'tdgpfa';
xDim = 8; 
kernSD = 30; % 30 (default) - for two-stage methods
segLength = 20; % 20 (default) 
binWidth = 20; % 20 (default)
result = neuralTraj(runIdx, dat, 'method', method, 'xDim', xDim,... 
                    'kernSDList', kernSD, 'segLength', segLength,...
                    'binWidth', binWidth);

plotDelayMatrix(result_tdgpfa.estParams, result_tdgpfa.binWidth);

plotTDGPFAlatentsVsTime(result_tdgpfa.seqTrain, result_tdgpfa.estParams, ...
    result_tdgpfa.binWidth);
                

if contains(fn,'JEB7')
    redTrials = find(obj.bp.L&obj.bp.hit&~obj.bp.early&(obj.bp.autowater.nums'==2)); % JEB7
else
    redTrials = find(obj.bp.L&obj.bp.hit&~obj.bp.early&~obj.bp.autowater); % EKH2
end
nPlotMax = 100;
plotEachDimVsTime(seqTrain, 'xorth', result.binWidth,'nPlotMax',nPlotMax,...
                  'redTrials',redTrials);
nPlotMax3D = 20;
plot3D(seqTrain, 'xorth', 'dimsToPlot', [3 4 5], 'nPlotMax',nPlotMax3D,...
                  'redTrials',redTrials);

'a'


