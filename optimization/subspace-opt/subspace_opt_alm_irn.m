%%subspace optimization code from https://github.com/jcykao/subspace-opt 
% Finds orthogonal subspaces for movement and preparatory neural activity
% calculates variance explained for each epoch and subspace
% plots projections of neural activity per condition onto subspaces
clear,clc,close all
% mandatory import of manopt functions
run('manopt/importmanopt'); % I turned off 'save path' question in importmanopt.m

% load data
data_pth = '/Users/Munib/Documents/Economo-Lab/subspace-id/data/elsayed/';
whos_data = 'erin';

data_file = 'EEL62020-03-01_EEL72020-03-02_processed_softnorm.mat'; % prep subspace looks ok
% data_file = 'EEL62020-03-01_EEL72020-03-02_processed_zscore.mat'; % prep subspace looks bad

load(fullfile(data_pth,whos_data,data_file));

currpath = pwd;
% this is where the optimization functions reside
addpath([currpath '/optFunctions']);

% path to store figures
fig_pth = '/Users/Munib/Documents/Economo-Lab/subspace-id/figs/elsayed/';
sav = 1; % save figs this run?


%% TODO
% sample equal populations from alm and irn
% identify subspaces
% repeat
% distribution of variance explained


%% Optimization 1, Orth Subspace for IRN and ALM jointly. Look at neuron contributions
d_Move = 4;
d_Prep = 4;

alpha = 0;  % regularization hyperparam (+ve->discourage sparity, -ve->encourage sparsity)

numIters = 100; % number of bootstrap iterations

% number of cells to sample each iteration
sampSize = round(0.9 * min([length(data.alm_idx),length(data.irn_idx)]));

move_on_move = zeros(1,numIters);
prep_on_move = zeros(1,numIters);
prep_on_prep = zeros(1,numIters);
move_on_prep = zeros(1,numIters);
alm_move_contrib = zeros(numIters,sampSize);
alm_prep_contrib = zeros(numIters,sampSize);
irn_move_contrib = zeros(numIters,sampSize);
irn_prep_contrib = zeros(numIters,sampSize);
for sub_iter = 1:numIters
    % sample from irn and alm
    almSampIdx = randi([data.alm_idx(1),data.alm_idx(end)],[1,sampSize]);
    irnSampIdx = randi([data.irn_idx(1),data.irn_idx(end)],[1,sampSize]);
    sampIdx = [almSampIdx,irnSampIdx];
    
    % % identify subspaces using sample
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
    % alm neuron contribution
    [alm_move_contrib(sub_iter,:),alm_prep_contrib(sub_iter,:)] = ...
                      subspace_contribution(Q(1:sampSize,:),P1,P2,data.mean_fr(almSampIdx));
    % irn neuron contribution 
    [irn_move_contrib(sub_iter,:),irn_prep_contrib(sub_iter,:)] = ...
                      subspace_contribution(Q((sampSize+1):sampSize*2,:),P1,P2,data.mean_fr(irnSampIdx));
                  
    % % identify subspaces using alm and irn separately
    % then project and cross-project
    % alm
    [Q_alm,~,P,~,~] = orthogonal_subspaces(data.Cmove(almSampIdx,almSampIdx),d_Move,...
                                           data.Cprep(almSampIdx,almSampIdx),d_Prep,alpha);
    P1_alm = P{1};
    P2_alm = P{2};

    % irn
    [Q_irn,~,P,~,~] = orthogonal_subspaces(data.Cmove(irnSampIdx,irnSampIdx),d_Move,...
                                           data.Cprep(irnSampIdx,irnSampIdx),d_Prep,alpha);
    P1_irn = P{1};
    P2_irn = P{2};
    
    % alm in alm subspaces
    % mm_aa = move in move for alm activity in alm subpsaces
    % mp_ai = move in prep for alm activity in irn subspaces
    [mm_aa(sub_iter),pm_aa(sub_iter),pp_aa(sub_iter),mp_aa(sub_iter)] =...
                      samp_var_explained(Q_alm,P1_alm,P2_alm,...
                      data.Cmove(almSampIdx,almSampIdx),...
                      data.Cprep(almSampIdx,almSampIdx),d_Move,d_Prep,dmax);

    % irn in irn subspaces
    [mm_ii(sub_iter),pm_ii(sub_iter),pp_ii(sub_iter),mp_ii(sub_iter)] =...
                      samp_var_explained(Q_irn,P1_irn,P2_irn,...
                      data.Cmove(irnSampIdx,irnSampIdx),...
                      data.Cprep(irnSampIdx,irnSampIdx),d_Move,d_Prep,dmax);
                  
    % alm in irn subspaces
    [mm_ai(sub_iter),pm_ai(sub_iter),pp_ai(sub_iter),mp_ai(sub_iter)] =...
                      samp_var_explained(Q_irn,P1_irn,P2_irn,...
                      data.Cmove(almSampIdx,almSampIdx),...
                      data.Cprep(almSampIdx,almSampIdx),d_Move,d_Prep,dmax);

    % irn in alm subspaces
    [mm_ia(sub_iter),pm_ia(sub_iter),pp_ia(sub_iter),mp_ia(sub_iter)] =...
                      samp_var_explained(Q_alm,P1_alm,P2_alm,...
                      data.Cmove(irnSampIdx,irnSampIdx),...
                      data.Cprep(irnSampIdx,irnSampIdx),d_Move,d_Prep,dmax);
    
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
fig_name = fullfile(fig_pth,whos_data,'joint_var_explained');
if sav; saveas(gcf,fig_name,'png'); end

% plot distribution of firing rates for each population
numBins = 10;
histAlpha = 0.3;

figure
histogram(data.mean_fr(data.alm_idx),numBins,'FaceAlpha',histAlpha,'Normalization','pdf'); hold on
histogram(data.mean_fr(data.irn_idx),numBins,'FaceAlpha',histAlpha,'Normalization','pdf'); hold on
xlabel('Mean Firing Rate')
title('Distribution of mean firing rate by group')
legend('ALM','IRN');
fig_name = fullfile(fig_pth,whos_data,'fr_dist');
if sav; saveas(gcf,fig_name,'png'); end

% plot histogram of contributions
numBins = round(sampSize / 4);
histAlpha = 0.3;

figure
subplot(1,2,1)
histogram(alm_move_contrib(:),numBins,'FaceAlpha',histAlpha); hold on
histogram(alm_prep_contrib(:),numBins,'FaceAlpha',histAlpha); hold off
xlabel('Contribution')
title('ALM Cells')
legend('Move','Prep')
subplot(1,2,2)
histogram(irn_move_contrib(:),numBins,'FaceAlpha',histAlpha); hold on
histogram(irn_prep_contrib(:),numBins,'FaceAlpha',histAlpha); hold off
xlabel('Contribution')
title('IRN Cells')
sgtitle('Single Neuron Subspace Contributions')
fig_name = fullfile(fig_pth,whos_data,'contribution_dist');
if sav; saveas(gcf,fig_name,'png'); end

% plot median contribution for each population to each subspace
figure
subplot(1,2,1)
bar([median(alm_move_contrib(:)),median(alm_prep_contrib(:))]);
grid on; ax = gca(); ax.XTickLabel = {'Move','Prep'};
a = get(gca,'XTickLabel'); set(gca,'XTickLabel',a,'fontsize',18)
xtickangle(45)
ylabel('Median Contribution');
title('ALM cells')
subplot(1,2,2)
bar([median(irn_move_contrib(:)),median(irn_prep_contrib(:))]);
grid on; ax = gca(); ax.XTickLabel = {'Move','Prep'};
a = get(gca,'XTickLabel'); set(gca,'XTickLabel',a,'fontsize',18)
xtickangle(45)
xlabel('Subspace');
title('IRN cells')
fig_name = fullfile(fig_pth,whos_data,'median_contribution');
if sav; saveas(gcf,fig_name,'png'); end

% variance explained (projections and cross-projections
numBins = round(numIters / 4);
histAlpha = 0.3;

figure
subplot(2,2,1)
histogram(mm_aa,numBins,'FaceAlpha',histAlpha); hold on
histogram(pp_aa,numBins,'FaceAlpha',histAlpha); hold on
histogram(pm_aa,numBins,'FaceAlpha',histAlpha); hold on
histogram(mp_aa,numBins,'FaceAlpha',histAlpha); hold off
legend('Move in Move','Prep in Prep','Prep in Move','Move in Prep')
title('ALM activity in ALM subspaces')
subplot(2,2,2)
histogram(mm_ii,numBins,'FaceAlpha',histAlpha); hold on
histogram(pp_ii,numBins,'FaceAlpha',histAlpha); hold on
histogram(pm_ii,numBins,'FaceAlpha',histAlpha); hold on
histogram(mp_ii,numBins,'FaceAlpha',histAlpha); hold off
legend('Move in Move','Prep in Prep','Prep in Move','Move in Prep')
title('IRN activity in IRN subspaces')
subplot(2,2,3)
histogram(mm_ai,numBins,'FaceAlpha',histAlpha); hold on
histogram(pp_ai,numBins,'FaceAlpha',histAlpha); hold on
histogram(pm_ai,numBins,'FaceAlpha',histAlpha); hold on
histogram(mp_ai,numBins,'FaceAlpha',histAlpha); hold off
legend('Move in Move','Prep in Prep','Prep in Move','Move in Prep')
title('ALM activity in IRN subspaces')
subplot(2,2,4)
histogram(mm_ia,numBins,'FaceAlpha',histAlpha); hold on
histogram(pp_ia,numBins,'FaceAlpha',histAlpha); hold on
histogram(pm_ia,numBins,'FaceAlpha',histAlpha); hold on
histogram(mp_ia,numBins,'FaceAlpha',histAlpha); hold off
legend('Move in Move','Prep in Prep','Prep in Move','Move in Prep')
title('IRN activity in ALM subspaces')

fig_name = fullfile(fig_pth,whos_data,'var_explained_irn_alm');
if sav; saveas(gcf,fig_name,'png'); end


%% Helper Functions

function var_explained(Q,P1,P2,C1,C2,d1,d2,dmax)
    eigvals1 = eigs(C1, dmax, 'la'); 
    eigvals2 = eigs(C2, dmax, 'la');
    Move_on_Move = var_proj(Q*P1,C1,sum(eigvals1(1:d1))); % var explained of Move in Orth-Move subsapce
    Prep_on_Move = var_proj(Q*P1,C2,sum(eigvals2(1:d1))); % var explained of Prep in Orth-Move subsapce
    Prep_on_Prep = var_proj(Q*P2,C2,sum(eigvals2(1:d2)));
    Move_on_Prep = var_proj(Q*P2,C1,sum(eigvals1(1:d2)));

    figure();
    bar([Prep_on_Prep, Prep_on_Move,0,0,Move_on_Move, Move_on_Prep]);
    grid on;
    ax = gca();
    ax.XTickLabel = {'Prep in Prep','Prep in Move','','','Move in Move','Move in Prep'};
    a = get(gca,'XTickLabel'); set(gca,'XTickLabel',a,'fontsize',18)
    xtickangle(45)
    xlabel('Subspace projections');
    ylabel('Fraction of variance captured');
end % var_explained

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

function [samp_move_contrib,samp_prep_contrib] = subspace_contribution(Q,P1,P2,fr)
    % psth = (TC,N) for a given population
    % contrib = ||Q_i||^2 * FR_i, where is a given neuron
    % each row of Q and column of psth corresponds to a single neuron
    samp_move_contrib = zeros(1,size(Q,1));
    samp_prep_contrib = zeros(1,size(Q,1));
    for i = 1:size(Q,1)
        samp_move_contrib(i) = norm(Q(i,:)*P1) * fr(i);
        samp_prep_contrib(i) = norm(Q(i,:)*P2) * fr(i);
    end

end % subspace_contribution