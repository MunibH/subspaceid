%% subspace optimization 
% code from https://github.com/jcykao/subspace-opt 

% Data details
% Population = 1:1207
% PT Lower   = 1208:1276
% PT Upper   = 1277:1337
% 
% -since each of these populations has unique number of cells, need to
% sample equal number of cells from each and then identify subspaces. This
% makes it so that the dimension of the subspaces are compatible with each
% of the populations. 
% -since the number of Population cells well exceeds number of PT Lower and
% PT Upper cells, each sample will be separately sampled from each
% population.

clear,clc,close all
% mandatory import of manopt functions
run('manopt/importmanopt'); % I turned off 'save path' question in importmanopt.m
warning('off', 'manopt:getHessian:approx')

% load data
data_pth = '/Users/Munib/Documents/Economo-Lab/subspace-id/data/elsayed/';
whos_data = 'mike';

% data_file = 'alldat_tagged_processed.mat'; % 1337 cells, soft-norm
% data_file = 'alldat_tagged_processed_zscore.mat';  % 1337 cells, z-score

load(fullfile(data_pth,whos_data,data_file));

currpath = pwd;
% this is where the optimization functions reside
addpath([currpath '/optFunctions']);

% path to store figures
fig_pth = '/Users/Munib/Documents/Economo-Lab/subspace-id/figs/elsayed/';
sav = 0; % save figs this run?

%% Method 1, Orth Subspace, Bootstrap over subspace identification
d_Move = 4;
d_Prep = 2;

% bootstrap over all population and tagged cells to identify subspace
numIters = 30;
% number of cells to sample each subspace iteration
sampSize = round(0.9 * min([length(data.pop_idx),length(data.ptlow_idx),length(data.ptup_idx)]));

move_on_move = zeros(1,numIters);
prep_on_move = zeros(1,numIters);
prep_on_prep = zeros(1,numIters);
move_on_prep = zeros(1,numIters);
pop_move_contrib = zeros(numIters,sampSize);
pop_prep_contrib = zeros(numIters,sampSize);
ptlow_move_contrib = zeros(numIters,sampSize);
ptlow_prep_contrib = zeros(numIters,sampSize);
ptup_move_contrib = zeros(numIters,sampSize);
ptup_prep_contrib = zeros(numIters,sampSize);
for sub_iter = 1:numIters
    % % sample from tagged+untagged population w/replacement
    popSampIdx = randi([data.pop_idx(1),data.pop_idx(end)],[1,sampSize]);
    ptlowSampIdx = randi([data.ptlow_idx(1),data.ptlow_idx(end)],[1,sampSize]);
    ptupSampIdx = randi([data.ptup_idx(1),data.ptup_idx(end)],[1,sampSize]);
    sampIdx = [popSampIdx,ptlowSampIdx,ptupSampIdx];
    
    % % identify subspaces using sample
    alpha = 0; % regularization hyperparam (+ve->discourage sparity, -ve->encourage sparsity)
    [Q,~,P,~,~] = orthogonal_subspaces(data.Cmove(sampIdx,sampIdx),d_Move,...
                                       data.Cprep(sampIdx,sampIdx),d_Prep,alpha);
    P1 = P{1}; % Qmove = Q*P1;
    P2 = P{2}; % Qprep = Q*P2;
    dmax = max(d_Move,d_Prep);
    
    % % variance captured for entire sampled population
    [move_on_move(sub_iter),prep_on_move(sub_iter),prep_on_prep(sub_iter),move_on_prep(sub_iter)] =... 
                                      samp_var_explained(Q,P1,P2,data.Cmove(sampIdx,sampIdx),...
                                      data.Cprep(sampIdx,sampIdx),d_Move,d_Prep,dmax);
                                  
    % % calculate each neurons contribution to the subspace for each
    % population neuron contribution
    [pop_move_contrib(sub_iter,:),pop_prep_contrib(sub_iter,:)] = ...
                      neuron_contribution(Q(1:sampSize,:),P1,P2,data.mean_fr(popSampIdx));
    % pt low neuron contribution 
    [ptlow_move_contrib(sub_iter,:),ptlow_prep_contrib(sub_iter,:)] = ...
                      neuron_contribution(Q((sampSize+1):sampSize*2,:),P1,P2,data.mean_fr(ptlowSampIdx));
    % pt low neuron contribution 
    [ptup_move_contrib(sub_iter,:),ptup_prep_contrib(sub_iter,:)] = ...
                      neuron_contribution(Q((sampSize*2+1):sampSize*3,:),P1,P2,data.mean_fr(ptupSampIdx));
  
    % % variance captured for individual sampled populations
    % non-tagged cells
    [move_on_move_nontag(sub_iter),prep_on_move_nontag(sub_iter),prep_on_prep_nontag(sub_iter),move_on_prep_nontag(sub_iter)] =... 
                                      samp_var_explained(Q(1:sampSize,:),P1,P2,data.Cmove(popSampIdx,popSampIdx),...
                                      data.Cprep(popSampIdx,popSampIdx),d_Move,d_Prep,dmax);
    % pt low
    [move_on_move_ptlow(sub_iter),prep_on_move_ptlow(sub_iter),prep_on_prep_ptlow(sub_iter),move_on_prep_ptlow(sub_iter)] =... 
                                      samp_var_explained(Q((sampSize+1):sampSize*2,:),P1,P2,data.Cmove(ptlowSampIdx,ptlowSampIdx),...
                                      data.Cprep(ptlowSampIdx,ptlowSampIdx),d_Move,d_Prep,dmax);
    % pt up
    [move_on_move_ptup(sub_iter),prep_on_move_ptup(sub_iter),prep_on_prep_ptup(sub_iter),move_on_prep_ptup(sub_iter)] =... 
                                      samp_var_explained(Q((sampSize*2+1):sampSize*3,:),P1,P2,data.Cmove(ptupSampIdx,ptupSampIdx),...
                                      data.Cprep(ptupSampIdx,ptupSampIdx),d_Move,d_Prep,dmax);
                                  
    % calculate whole population contributions (nontagged,ptlow,ptup)
    % % calculate each neurons contribution to the subspace for each
    % population neuron contribution
    [all_nontag_move_contrib(sub_iter),all_nontag_prep_contrib(sub_iter)] = ...
                      population_contribution(Q(1:sampSize,:),P1,P2,data.mean_fr(popSampIdx));
    % pt low neuron contribution 
    [all_ptlow_move_contrib(sub_iter),all_ptlow_prep_contrib(sub_iter)] = ...
                      population_contribution(Q((sampSize+1):sampSize*2,:),P1,P2,data.mean_fr(ptlowSampIdx));
    % pt low neuron contribution 
    [all_ptup_move_contrib(sub_iter),all_ptup_prep_contrib(sub_iter)] = ...
                      population_contribution(Q((sampSize*2+1):sampSize*3,:),P1,P2,data.mean_fr(ptupSampIdx));
    
end

%% Plotting 

% plot variance explained for bootstrap sampling distribution
numBins = round(numIters / 4);
histAlpha = 0.3;

figure
subplot(2,1,1)
histogram(move_on_move,numBins,'FaceAlpha',histAlpha); hold on
histogram(prep_on_prep,numBins,'FaceAlpha',histAlpha); hold off
xlim([0.5,1])
legend('Move in Move','Prep in Prep')
subplot(2,1,2)
histogram(prep_on_move,numBins,'FaceColor',[0.9290 0.6940 0.1250],'FaceAlpha',histAlpha); hold on
histogram(move_on_prep,numBins,'FaceColor',[0.4940 0.1840 0.5560],'FaceAlpha',histAlpha); hold off
xlim([0,0.5])
xlabel('Variance captured')
legend('Prep in Move','Move in Prep')
sgtitle('Sampling Distribution - variance explained')
fig_name = fullfile(fig_pth,whos_data,'var_explained');
if sav; saveas(gcf,fig_name,'png'); end

% plot distribution of firing rates for each population
numBins = 10;
histAlpha = 0.3;

figure
histogram(data.mean_fr(data.pop_idx),numBins,'FaceAlpha',histAlpha,'Normalization','pdf'); hold on
histogram(data.mean_fr(data.ptlow_idx),numBins,'FaceAlpha',histAlpha,'Normalization','pdf'); hold on
histogram(data.mean_fr(data.ptup_idx),numBins,'FaceAlpha',histAlpha,'Normalization','pdf'); hold on
xlabel('Mean Firing Rate')
title('Distribution of mean firing rate by group')
legend('Non-tagged','PT Lower','PT Upper');
fig_name = fullfile(fig_pth,whos_data,'fr_dist');
if sav; saveas(gcf,fig_name,'png'); end

% plot histogram of single neuron contributions
numBins = round(sampSize / 4);
histAlpha = 0.3;

figure
subplot(1,3,1)
histogram(pop_move_contrib(:),numBins,'FaceAlpha',histAlpha); hold on
histogram(pop_prep_contrib(:),numBins,'FaceAlpha',histAlpha); hold off
xlabel('Contribution')
title('Non-tagged Cells')
legend('Move','Prep')
subplot(1,3,2)
histogram(ptlow_move_contrib(:),numBins,'FaceAlpha',histAlpha); hold on
histogram(ptlow_prep_contrib(:),numBins,'FaceAlpha',histAlpha); hold off
xlabel('Contribution')
title('PT Lower Cells')
subplot(1,3,3)
histogram(ptup_move_contrib(:),numBins,'FaceAlpha',histAlpha); hold on
histogram(ptup_prep_contrib(:),numBins,'FaceAlpha',histAlpha); hold off
xlabel('Contribution')
title('PT Upper Cells')
sgtitle('Single Neuron Subspace Contributions')
fig_name = fullfile(fig_pth,whos_data,'contribution_dist');
if sav; saveas(gcf,fig_name,'png'); end

% plot median contribution for each population to each subspace
figure
subplot(1,3,1)
bar([median(pop_move_contrib(:)),median(pop_prep_contrib(:))]);
grid on; ax = gca(); ax.XTickLabel = {'Move','Prep'};
a = get(gca,'XTickLabel'); set(gca,'XTickLabel',a,'fontsize',18)
xtickangle(45)
ylabel('Median Contribution');
title('Non-tagged cells')
subplot(1,3,2)
bar([median(ptlow_move_contrib(:)),median(ptlow_prep_contrib(:))]);
grid on; ax = gca(); ax.XTickLabel = {'Move','Prep'};
a = get(gca,'XTickLabel'); set(gca,'XTickLabel',a,'fontsize',18)
xtickangle(45)
xlabel('Subspace');
title('PT Lower cells')
subplot(1,3,3)
bar([median(ptup_move_contrib(:)),median(ptup_prep_contrib(:))]);
grid on; ax = gca(); ax.XTickLabel = {'Move','Prep'};
a = get(gca,'XTickLabel'); set(gca,'XTickLabel',a,'fontsize',18)
xtickangle(45)
title('PT Upper cells')
fig_name = fullfile(fig_pth,whos_data,'median_contribution');
if sav; saveas(gcf,fig_name,'png'); end

% plot variance explained for separate bootstrap sampling distribution
numBins = round(numIters / 4);
histAlpha = 0.3;

figure
subplot(1,3,1)
histogram(move_on_move_nontag,numBins,'FaceAlpha',histAlpha); hold on
histogram(prep_on_prep_nontag,numBins,'FaceAlpha',histAlpha); hold on
histogram(prep_on_move_nontag,numBins,'FaceAlpha',histAlpha); hold on
histogram(move_on_prep_nontag,numBins,'FaceAlpha',histAlpha); hold off
title('Non-tagged')
subplot(1,3,2)
histogram(move_on_move_ptlow,numBins,'FaceAlpha',histAlpha); hold on
histogram(prep_on_prep_ptlow,numBins,'FaceAlpha',histAlpha); hold on
histogram(prep_on_move_ptlow,numBins,'FaceAlpha',histAlpha); hold on
histogram(move_on_prep_ptlow,numBins,'FaceAlpha',histAlpha); hold off
xlabel('Variance captured'); title('PT Low')
legend('Move in Move','Prep in Prep','Prep in Move','Move in Prep')
subplot(1,3,3)
histogram(move_on_move_ptup,numBins,'FaceAlpha',histAlpha); hold on
histogram(prep_on_prep_ptup,numBins,'FaceAlpha',histAlpha); hold on
histogram(prep_on_move_ptup,numBins,'FaceAlpha',histAlpha); hold on
histogram(move_on_prep_ptup,numBins,'FaceAlpha',histAlpha); hold off
title('PT Up')
sgtitle('variance explained by cell type')

% plot histogram of cell type contributions
numBins = round(numIters / 4);
histAlpha = 0.3;

figure
subplot(1,3,1)
histogram(all_nontag_move_contrib,numBins,'FaceAlpha',histAlpha); hold on
histogram(all_nontag_prep_contrib,numBins,'FaceAlpha',histAlpha); hold off
xlabel('Contribution')
title('Non-tagged Cells')
legend('Move','Prep')
subplot(1,3,2)
histogram(all_ptlow_move_contrib,numBins,'FaceAlpha',histAlpha); hold on
histogram(all_ptlow_prep_contrib,numBins,'FaceAlpha',histAlpha); hold off
xlabel('Contribution')
title('PT Lower Cells')
subplot(1,3,3)
histogram(all_ptup_move_contrib,numBins,'FaceAlpha',histAlpha); hold on
histogram(all_ptup_prep_contrib,numBins,'FaceAlpha',histAlpha); hold off
xlabel('Contribution')
title('PT Upper Cells')
sgtitle('Cell Type Subspace Contributions')
fig_name = fullfile(fig_pth,whos_data,'contribution_dist');
if sav; saveas(gcf,fig_name,'png'); end


% %% Method 2, Orth Subspace, Bootstrap over variance explained
% d_Move = d_Move;
% d_Prep = d_Prep;
% 
% % bootstrap over all population and tagged cells to identify subspace
% numIters = 1000;
% 
% % identify subspaces with all cells (tagged + untagged)
% alpha = 0; % regularization hyperparam (+ve->discourage sparity, -ve->encourage sparsity)
% [Q,~,P,~,~] = orthogonal_subspaces(data.Cmove,d_Move,...
%                                    data.Cprep,d_Prep,alpha);
% P1 = P{1}; % Qmove = Q*P1;
% P2 = P{2}; % Qprep = Q*P2;
% dmax = max(d_Move,d_Prep);
% 
% % sample size
% sampSize = round(0.9 * min([length(data.pop_idx),length(data.ptlow_idx),length(data.ptup_idx)]));
% % random cells
% cells_to_use = [1,size(data.Cprep,1)];
% bootstrap_var_explained(Q,P1,P2,data.Cmove,data.Cprep,d_Move,d_Prep,dmax,numIters,sampSize,cells_to_use)
% title('random')
% fig_name = fullfile(fig_pth,whos_data,'random');
% if sav; saveas(gcf,fig_name,'png'); end
% % non tagged cells
% cells_to_use = [data.pop_idx(1),data.pop_idx(end)];
% bootstrap_var_explained(Q,P1,P2,data.Cmove,data.Cprep,d_Move,d_Prep,dmax,numIters,sampSize,cells_to_use)
% title('non-tagged')
% fig_name = fullfile(fig_pth,whos_data,'nontagged');
% if sav; saveas(gcf,fig_name,'png'); end
% % pt lower cells
% cells_to_use = [data.ptlow_idx(1),data.ptlow_idx(end)];
% bootstrap_var_explained(Q,P1,P2,data.Cmove,data.Cprep,d_Move,d_Prep,dmax,numIters,sampSize,cells_to_use)
% title('PT Lower')
% fig_name = fullfile(fig_pth,whos_data,'ptlow');
% if sav; saveas(gcf,fig_name,'png'); end
% % pt upper cells
% cells_to_use = [data.ptup_idx(1),data.ptup_idx(end)];
% bootstrap_var_explained(Q,P1,P2,data.Cmove,data.Cprep,d_Move,d_Prep,dmax,numIters,sampSize,cells_to_use)
% title('PT Upper')
% fig_name = fullfile(fig_pth,whos_data,'ptup');
% if sav; saveas(gcf,fig_name,'png'); end

%% Helper Functions

function bootstrap_var_explained(Q,P1,P2,Cmove,Cprep,d_Move,d_Prep,dmax,numIters,nCells,cell_ix)
    move_on_move_dist = zeros(1,numIters);
    prep_on_prep_dist = zeros(1,numIters);
    prep_on_move_dist = zeros(1,numIters);
    move_on_prep_dist = zeros(1,numIters);
    for i = 1:numIters
        ix = randi(cell_ix,[1,nCells]);
        eg1 = eigs(Cmove(ix,ix), dmax, 'la'); 
        eg2 = eigs(Cprep(ix,ix), dmax, 'la');
        Move_on_Move = var_proj(Q(ix,:)*P1,Cmove(ix,ix),sum(eg1(1:d_Move))); % var explained of Move in Orth-Move subsapce
        Prep_on_Move = var_proj(Q(ix,:)*P1,Cprep(ix,ix),sum(eg2(1:d_Move))); % var explained of Prep in Orth-Move subsapce
        Prep_on_Prep = var_proj(Q(ix,:)*P2,Cprep(ix,ix),sum(eg2(1:d_Prep)));
        Move_on_Prep = var_proj(Q(ix,:)*P2,Cmove(ix,ix),sum(eg1(1:d_Prep)));
        move_on_move_dist(i) = Move_on_Move;
        prep_on_prep_dist(i) = Prep_on_Prep;
        prep_on_move_dist(i) = Prep_on_Move;
        move_on_prep_dist(i) = Move_on_Prep;
    end
    numBins = 20;
    figure
    histogram(move_on_move_dist,numBins,'FaceAlpha',0.3); xlim([0,0.2]);
    hold on
    histogram(prep_on_prep_dist,numBins,'FaceAlpha',0.3); xlim([0,0.2]);
    hold on
    histogram(prep_on_move_dist,numBins,'FaceAlpha',0.3); xlim([0,0.2]);
    hold on
    histogram(move_on_prep_dist,numBins,'FaceAlpha',0.3); xlim([0,0.2]);
    xlabel('Variance captured')
    legend('Move in Move','Prep in Prep')
end % bootstrap_var_explained

function [move_on_move,prep_on_move,prep_on_prep,move_on_prep] = samp_var_explained(Q,P1,P2,C1,C2,d1,d2,dmax)
    eigvals1 = eigs(C1, dmax, 'la'); 
    eigvals2 = eigs(C2, dmax, 'la');
    move_on_move = var_proj(Q*P1,C1,sum(eigvals1(1:d1))); % var explained of Move in Orth-Move subsapce
    prep_on_move = var_proj(Q*P1,C2,sum(eigvals2(1:d1))); % var explained of Prep in Orth-Move subsapce
    prep_on_prep = var_proj(Q*P2,C2,sum(eigvals2(1:d2)));
    move_on_prep = var_proj(Q*P2,C1,sum(eigvals1(1:d2)));
end  % samp_var_explained

function [samp_move_contrib,samp_prep_contrib] = neuron_contribution(Q,P1,P2,fr)
    % psth = (TC,N) for a given population
    % contrib = ||Q_i||^2 * FR_i, where is a given neuron
    % each row of Q and column of psth corresponds to a single neuron
    samp_move_contrib = zeros(1,size(Q,1));
    samp_prep_contrib = zeros(1,size(Q,1));
    for i = 1:size(Q,1)
        samp_move_contrib(i) = norm(Q(i,:)*P1) * fr(i);
        samp_prep_contrib(i) = norm(Q(i,:)*P2) * fr(i);
    end

end % neuron_contribution

function [samp_move_contrib,samp_prep_contrib] = population_contribution(Q,P1,P2,fr)
    % psth = (TC,N) for a given population
    % contrib = ||Q_i||^2 * FR_i, where is a given neuron
    % each row of Q and column of psth corresponds to a single neuron
%     samp_move_contrib = norm(Q*P1,'fro') * mean(fr);
%     samp_prep_contrib = norm(Q*P2,'fro') * mean(fr);
    samp_move_contrib = zeros(1,size(Q,1));
    samp_prep_contrib = zeros(1,size(Q,1));
    for i = 1:size(Q,1)
        samp_move_contrib(i) = norm(Q(i,:)*P1) * fr(i);
        samp_prep_contrib(i) = norm(Q(i,:)*P2) * fr(i);
    end
    samp_move_contrib = sum(samp_move_contrib);
    samp_prep_contrib = sum(samp_prep_contrib);
    

end % population_contribution











