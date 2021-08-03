%%subspace optimization code from https://github.com/jcykao/subspace-opt 
% Finds orthogonal subspaces for movement and preparatory neural activity
% calculates variance explained for each epoch and subspace
% plots projections of neural activity per condition onto subspaces
clear,clc,close all
% mandatory import of manopt functions
run('manopt/importmanopt'); % I turned off 'save path' question in importmanopt.m
warning('off', 'manopt:getHessian:approx')

% load data
data_pth = '/Users/Munib/Documents/Economo-Lab/subspace-id/data/elsayed/';
whos_data = 'mike';

% data_file = 'alldat_processed.mat'; % 1207 cells
% data_file = 'alldat_tagged_processed.mat'; % 1337 cells
% data_file = 'alldat_tagged_processed_zscore.mat';  % 1337 cells, z-score
% data_file = 'alldat_tagged_processed_mcwithin.mat'; % 1337 cells, mean-centered within conditions
% data_file = 'alldat_tagged_processed_int.mat'; % 1337 cells, soft-norm, contains psth_int that is the regular mean-centered psths
data_file = 'alldat_tagged_for_grant.mat'; % 1337 cells, extra smooth
load(fullfile(data_pth,whos_data,data_file));

currpath = pwd;
% this is where the optimization functions reside
addpath([currpath '/optFunctions']);

% path to store figures
fig_pth = '/Users/Munib/Documents/Economo-Lab/subspace-id/figs/elsayed/';
sav = 0; % save figs this run?

%% Optimization 1, Orth Subspace
d_Move = 4;
d_Prep = 4;
dmax = max(d_Move,d_Prep);

[Q_null,~,~,~] = null_subspace(data.Cprep,d_Prep);
% [Q_potent,~,~,~] = potent_subspace_constrained(Q_null,data.Cmove,d_Prep,data.Cprep,d_Prep);
[Q_potent,~,~,~] = potent_subspace(Q_null,data.Cmove,d_Prep,data.Cprep,d_Prep);

P1 = [eye(d_Move); zeros(d_Prep,d_Move)];
P2 = [zeros(d_Move, d_Prep); eye(d_Prep)];

% % use gram-schmidt process to orthogonalize Q_potent and Q_null
% this process works because Q_null doesn't change and
% the columns of Q_potent are orthogonalized to the columns of Q_null
% and to the other columns of Q_potent
Q = gschmidt([Q_null , Q_potent]);
tempQ = Q;
if dmax == 4
    tempQ(:,1:4) = tempQ(:,5:8);
    tempQ(:,5:8) = Q(:,1:4);
elseif dmax == 2
    tempQ(:,1:2) = tempQ(:,3:4);
    tempQ(:,3:4) = Q(:,1:2);
end
Q = tempQ;

a = [1;0;0]; b = [0;0;1];
test = gschmidt([a,b]);

%% variance explained 

% for each condition in each subspace
var_explained(Q,P2,P1,data.Cmove,data.Cprep,d_Move,d_Prep,dmax);
title('Variance captured');
fig_name = fullfile(fig_pth,whos_data,[data_file '_var_explained']);
if sav; saveas(gcf,fig_name,'png'); end


%% projections for each condition in each subspace dimension
Q_move = Q*P2;
Q_prep = Q*P1;

full_left_move = data.full_psth{1}*Q_move;
full_right_move = data.full_psth{2}*Q_move;
full_left_prep = data.full_psth{1}*Q_prep;
full_right_prep = data.full_psth{2}*Q_prep;

% plot full trial projections
% move subspace
x_lim = [data.full_time(1),data.full_time(end)];
y_lim(1) = min( [min(min(full_left_move)), min(min(full_right_move))] );
y_lim(2) = max( [max(max(full_left_move)), max(max(full_right_move))] );
plot_projections(data.full_time,full_left_move,full_right_move,x_lim,y_lim,d_Move,'Move')
fig_name = fullfile(fig_pth,whos_data,[data_file '_move_proj']);
if sav; saveas(gcf,fig_name,'png'); end

% prep subspace
x_lim = [data.full_time(1),data.full_time(end)];
y_lim(1) = min( [min(min(full_left_prep)), min(min(full_right_prep))] );
y_lim(2) = max( [max(max(full_left_prep)), max(max(full_right_prep))] );
plot_projections(data.full_time,full_left_prep,full_right_prep,x_lim,y_lim,d_Prep,'Prep')
fig_name = fullfile(fig_pth,whos_data,[data_file '_prep_proj']);
if sav; saveas(gcf,fig_name,'png'); end

%% state trajectories

% plot full trial in state space
dt = mode(diff(data.full_time));
[~,t_go_idx] = min(abs(data.full_time - 0));
prep_sidx = floor(t_go_idx - (1.15 / dt)); % start of prep
prep_eidx = floor(t_go_idx - (0.15 / dt)); % end of prep
move_sidx = floor(t_go_idx + 1);
move_eidx = floor(t_go_idx + (1.0 / dt));

% colors {early,prep,move,late,prep start,go cue,move end}
right_cols = {1/255*[137, 207, 240],1/255*[102, 153, 204],1/255*[25, 25, 112],1/255*[0, 0, 139],1/255*[0, 255, 209],1/255*[31, 253, 68],1/255*[1, 188, 32]}; % blue,cool colors
left_cols = {1/255*[234, 60, 83],1/255*[180, 55, 87],1/255*[191, 10, 48],1/255*[150, 0, 24],1/255*[255, 213, 0],1/255*[255, 152, 0],1/255*[255, 67, 0]}; % red,warm colors

plot_full_state_space(full_left_move,full_left_prep,full_right_move,full_right_prep,...
                      t_go_idx,prep_sidx,prep_eidx,move_sidx,move_eidx,right_cols,left_cols);

% my_animate(gcf,fig_pth,whos_data);

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
    ax.XTickLabel = {'Prep in Null','Prep in Potent','','','Move in Potent','Move in Null'};
    a = get(gca,'XTickLabel'); set(gca,'XTickLabel',a,'fontsize',18)
    xtickangle(45)
    xlabel('Subspace projections');
    ylabel('Fraction of variance captured');
end % var_explained

function bootstrap_var_explained(Q,P1,P2,Cmove,Cprep,d_Move,d_Prep,dmax,numIters,nCells,cell_ix)
    move_on_move_dist = zeros(1,numIters);
    prep_on_prep_dist = zeros(1,numIters);
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
    figure
    histogram(move_on_move_dist,40,'FaceAlpha',0.3)
    hold on
    histogram(prep_on_prep_dist,40,'FaceAlpha',0.3)
    xlabel('Variance captured')
    legend('Move in Move','Prep in Prep')
end % bootstrap_var_explained

function plot_projections(time,left,right,x_lim,y_lim,d,title_str)
    figure
    for i = 1:d
        subplot(d,1,i)
%         figure
        plot(time,left(:,i),'r','LineWidth',2); hold on
        plot(time,right(:,i),'b','LineWidth',2); xlim(x_lim); ylim(y_lim); hold on
        set(gca,'yticklabel',[]);
        set(gca,'fontsize',18); xlabel('time relative to go cue');
        if i~=d; set(gca,'xticklabel',[]); 
        else; set(gca,'fontsize',18); xlabel('time relative to go cue'); end
        fill([-1.15,-0.05,-0.05,-1.15],[y_lim(1),y_lim(1),y_lim(2),y_lim(2)],'r','FaceAlpha',0.05,'EdgeColor','none'); % fill prep epoch
        hold on 
        fill([0.05,1.15,1.15,0.05],[y_lim(1),y_lim(1),y_lim(2),y_lim(2)],'g','FaceAlpha',0.05,'EdgeColor','none'); % fill move epoch
        hold off
        sgtitle([title_str ' Dim ' num2str(i)])
    end
end

function plot_full_state_space(full_left_move,full_left_prep,full_right_move,full_right_prep,...
                      t_go_idx,prep_sidx,prep_eidx,move_sidx,move_eidx,right_cols,left_cols)
    figure
    ax = gca();
    legend('AutoUpdate','off')
    p1 = plot3(full_left_move(1:prep_sidx,1),full_left_move(1:prep_sidx,2),full_left_prep(1:prep_sidx,1),'Color',left_cols{1},'LineWidth',2,'DisplayName','left early');
    hold on
    p2 = plot3(full_left_move(prep_sidx:t_go_idx,1),full_left_move(prep_sidx:t_go_idx,2),full_left_prep(prep_sidx:t_go_idx,1),'Color',left_cols{2},'LineWidth',2,'DisplayName','left prep');
    hold on
    p3 = plot3(full_left_move(t_go_idx:move_eidx,1),full_left_move(t_go_idx:move_eidx,2),full_left_prep(t_go_idx:move_eidx,1),'Color',left_cols{3},'LineWidth',2,'DisplayName','left move');
    hold on
    p4 = plot3(full_left_move(move_eidx:end,1),full_left_move(move_eidx:end,2),full_left_prep(move_eidx:end,1),'Color',left_cols{4},'LineWidth',2,'DisplayName','left late');
    hold on
    p5 = plot3(full_left_move(prep_sidx,1),full_left_move(prep_sidx,2),full_left_prep(prep_sidx,1),'o','MarkerSize',12,'Color','k','MarkerFaceColor',left_cols{5},'DisplayName','prep start');
    hold on
    p6 = plot3(full_left_move(t_go_idx,1),full_left_move(t_go_idx,2),full_left_prep(t_go_idx,1),'o','MarkerSize',12,'Color','k','MarkerFaceColor',left_cols{6},'DisplayName','go cue');
    hold on
    p7 = plot3(full_left_move(move_eidx,1),full_left_move(move_eidx,2),full_left_prep(move_eidx,1),'o','MarkerSize',12,'Color','k','MarkerFaceColor',left_cols{7},'DisplayName','move end');
    hold on
    p8 = plot3(full_right_move(1:prep_sidx,1),full_right_move(1:prep_sidx,2),full_right_prep(1:prep_sidx,1),'Color',right_cols{1},'LineWidth',2,'DisplayName','right early');
    hold on
    p9 = plot3(full_right_move(prep_sidx:t_go_idx,1),full_right_move(prep_sidx:t_go_idx,2),full_right_prep(prep_sidx:t_go_idx,1),'Color',right_cols{2},'LineWidth',2,'DisplayName','right prep');
    hold on
    p10 = plot3(full_right_move(t_go_idx:move_eidx,1),full_right_move(t_go_idx:move_eidx,2),full_right_prep(t_go_idx:move_eidx,1),'Color',right_cols{3},'LineWidth',2,'DisplayName','right move');
    hold on
    p11 = plot3(full_right_move(move_eidx:end,1),full_right_move(move_eidx:end,2),full_right_prep(move_eidx:end,1),'Color',right_cols{4},'LineWidth',2,'DisplayName','right late');
    hold on
    p12 = plot3(full_right_move(prep_sidx,1),full_right_move(prep_sidx,2),full_right_prep(prep_sidx,1),'o','MarkerSize',12,'Color','k','MarkerFaceColor',right_cols{5},'DisplayName','prep start');
    hold on
    p13 = plot3(full_right_move(t_go_idx,1),full_right_move(t_go_idx,2),full_right_prep(t_go_idx,1),'o','MarkerSize',12,'Color','k','MarkerFaceColor',right_cols{6},'DisplayName','go cue');
    hold on
    p14 = plot3(full_right_move(move_eidx,1),full_right_move(move_eidx,2),full_right_prep(move_eidx,1),'o','MarkerSize',12,'Color','k','MarkerFaceColor',right_cols{7},'DisplayName','move end');

    grid on
    legend();
    xlabel('Move Dim 1')
    ylabel('Move Dim 2')
    zlabel('Prep Dim 1')
end % plot_full_state_space

function [Q,R] = gschmidt(V)
% Input: V is an m by n matrix of full rank m<=n
% Output: an m-by-n upper triangular matrix R
% and an m-by-m unitary matrix Q so that A = Q*R.
    [m,n]=size(V);
    R=zeros(n);
    R(1,1)=norm(V(:,1));
    Q(:,1)=V(:,1)/R(1,1);
    for k=2:n
        R(1:k-1,k)=Q(:,1:k-1)'*V(:,k);
        Q(:,k)=V(:,k)-Q(:,1:k-1)*R(1:k-1,k);
        R(k,k)=norm(Q(:,k));
        Q(:,k)=Q(:,k)/R(k,k);
    end
end % gschmidt

