%%subspace optimization code from https://github.com/jcykao/subspace-opt 
% Finds orthogonal subspaces for movement and preparatory neural activity
% calculates variance explained for each epoch and subspace
% plots projections of neural activity per condition onto subspaces
close all
% mandatory import of manopt functions
run('manopt/importmanopt'); % I turned off 'save path' question in importmanopt.m

% load the psths for each condition and covariance matrices from both epochs. 
data_pth = '/Users/Munib/Documents/Economo-Lab/subspace-id/data/elsayed/';

% whos_data = 'mike';
whos_data = 'erin';

% loads N_left,N_right,Cprep,Cmove,full_psth,full_time,ptlow_ix,ptup_ix
data_file = 'psths_cov.mat'; % all cells (including tagged)
load(fullfile(data_pth,whos_data,data_file));

currpath = pwd;
% this is where the optimization functions reside
addpath([currpath '/optFunctions']);

% path to store figures
fig_pth = '/Users/Munib/Documents/Economo-Lab/subspace-id/figs/elsayed/';
sav = 0; % save figs this run?

%% Optimization 1, Orth Subspace
d_Move = 4;
d_Prep = 2;

alpha = 0; % regularization hyperparam (+ve->discourage sparity, -ve->encourage sparsity)
[Q, ~, info, options] = orthogonal_subspaces(Cmove,d_Move,Cprep,d_Prep,alpha);
P1 = [eye(d_Move); zeros(d_Prep,d_Move)];
P2 = [zeros(d_Move, d_Prep); eye(d_Prep)];
dmax = max(d_Move,d_Prep);

%% variance explained 

% for each condition in each subspace
var_explained(Q,P1,P2,Cmove,Cprep,d_Move,d_Prep,dmax);
title('Variance captured for orthogonal subspace hypothesis');
fig_name = fullfile(fig_pth,whos_data,'tagged_var_explained_full');
if sav; saveas(gcf,fig_name,'png'); end

%% bootstrap sampling distribution of variance explained 

bootstrap = 0;
if bootstrap
    % of random cells (including tagged)
    numIters = 1000;
    nCells = 61; % smallest group size (pt_up cells)
    cell_to_use = size(Cmove,1); % all cells
    bootstrap_var_explained(Q,P1,P2,Cmove,Cprep,d_Move,d_Prep,dmax,numIters,nCells,cell_to_use);
    title(['Bootstrap distribution of var explained for random cells, n=' num2str(nCells)])
    fig_name = fullfile(fig_pth,whos_data,'bootstrap_var_explained_all');
    if sav; saveas(gcf,fig_name,'png'); end
end

%% projections for each condition in each subspace dimension
Q_move = Q*P1;
Q_prep = Q*P2;

left_move = (N_left*Q_move);
right_move = (N_right*Q_move);
left_prep = (N_left*Q_prep);
right_prep = (N_right*Q_prep);

full_left_move = (full_psth{1}*Q_move);
full_right_move = full_psth{2}*Q_move;
full_left_prep = full_psth{1}*Q_prep;
full_right_prep = full_psth{2}*Q_prep;

% plot full trial projections
% move subspace
x_lim = [min(full_time),full_time(end)];
y_lim = [-30,30];
plot_projections(full_time,full_left_move,full_right_move,x_lim,y_lim,d_Move,'Move')
fig_name = fullfile(fig_pth,whos_data,'full_move_projection');
if sav; saveas(gcf,fig_name,'png'); end

% prep subspace
x_lim = [min(full_time),full_time(end)];
y_lim = [-20,20];
plot_projections(full_time,full_left_prep,full_right_prep,x_lim,y_lim,d_Prep,'Prep')
fig_name = fullfile(fig_pth,whos_data,'full_prep_projection');
if sav; saveas(gcf,fig_name,'png'); end

%% state trajectories

% plot epochs in state space
prep_t = 1:(size(left_move,1)/2);
move_t = (prep_t(end)):size(left_move,1);
figure
plot3(left_move(prep_t,1),left_move(prep_t,2),left_prep(prep_t,1),'r-','LineWidth',0.5)
hold on
plot3(left_move(move_t,1),left_move(move_t,2),left_prep(move_t,1),'r-','LineWidth',3)
hold on
plot3(left_move(1,1),left_move(1,2),left_prep(1,1),'r.','MarkerSize',40)
hold on
plot3(left_move(end,1),left_move(end,2),left_prep(end,1),'r>','MarkerSize',20)
hold on
plot3(right_move(prep_t,1),right_move(prep_t,2),right_prep(prep_t,1),'b-','LineWidth',1);
hold on
plot3(right_move(move_t,1),right_move(move_t,2),right_prep(move_t,1),'b-','LineWidth',3);
hold on
plot3(right_move(1,1),right_move(1,2),right_prep(1,1),'b.','MarkerSize',40);
hold on
plot3(right_move(end,1),right_move(end,2),right_prep(end,1),'b>','MarkerSize',20);
grid on
legend('left prep','left move','prep start','move end','right prep','right move','prep start','move end')
xlabel('Move Dim 1')
ylabel('Move Dim 2')
zlabel('Prep Dim 1')

% plot full trial in state space
dt = mode(diff(full_time));
[~,t_go_idx] = min(abs(full_time - 0));
prep_sidx = floor(t_go_idx - (1.15 / dt)); % start of prep
prep_eidx = floor(t_go_idx - (0.15 / dt)); % end of prep
move_sidx = floor(t_go_idx + 1);
move_eidx = floor(t_go_idx + (1.0 / dt));

% colors {early,prep,move,late,prep start,go cue,move end}
right_cols = {1/255*[137, 207, 240],1/255*[102, 153, 204],1/255*[25, 25, 112],1/255*[0, 0, 139],1/255*[0, 255, 209],1/255*[31, 253, 68],1/255*[1, 188, 32]}; % blue,cool colors
left_cols = {1/255*[234, 60, 83],1/255*[180, 55, 87],1/255*[191, 10, 48],1/255*[150, 0, 24],1/255*[255, 213, 0],1/255*[255, 152, 0],1/255*[255, 67, 0]}; % red,warm colors

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

% my_animate(gcf,fig_pth,whos_data);


%% Linear Decoder (Prep explaining Move)

% prep_t = 101:201; % last half second of prep epoch
% move_t = 301:401; % last half second of move epoch, each epoch is 1 sec
% 
% Xprep = [left_prep(prep_t,:) right_prep(prep_t,:)];
% Xmove = [left_move(move_t,:) right_move(move_t,:)];
% 
% % least squares
% W = (Xprep'*Xprep)\(Xprep'*Xmove);
% r2 = 1 - (norm(Xmove - (Xprep*W),'fro')^2 / norm(Xmove,'fro')^2);
% 
% % plot
% figure();
% bar([0.9976,0.9504]) % from separate runs of mike and erin's data
% grid on;
% ax = gca();
% ax.XTickLabel = {'Mike','Erin'};
% a = get(gca,'XTickLabel'); set(gca,'XTickLabel',a,'fontsize',18)
% ylabel('R^2');
% title('Linear Decoding - explain move with prep subspace');
% fig_name = fullfile(fig_pth,whos_data,'linear_decoder');
% % if sav; saveas(gcf,fig_name,'png'); end

%% plot covariance matrices and psths

% % psths
% figure
% subplot(1,2,1);
% plot(time,N_left); xlim([-1,1]); ylim([-3,3]); title('Left')
% a = get(gca,'XTickLabel'); set(gca,'xticklabel',a,'fontsize',18);
% xlabel('time relative to go cue')
% xlabel('time relative to go cue')
% subplot(1,2,2);
% plot(time,N_right); xlim([-1,1]); ylim([-3,3]); title('Right')
% xlabel('time relative to go cue')
% sgtitle('psth - centered')
% a = get(gca,'XTickLabel'); set(gca,'XTickLabel',a,'fontsize',18)
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
% fig_name = fullfile(fig_pth,whos_data,'psth');
% % if sav; saveas(gcf,fig_name,'png'); end

% % psth heatmap
% % sort the psth's by mean activity in move epoch
% move_t = 201:402;
% [~,sortIdx] = sort(mean(N_left(move_t,:)));
% N_left_sorted = N_left(:,sortIdx);
% [~,sortIdx] = sort(mean(N_right(move_t,:)));
% N_right_sorted = N_right(:,sortIdx);
% 
% cutout = 0.05;
% figure
% colormap winter
% subplot(1,2,1);
% plot_cov(N_left_sorted',cutout); title('Left')
% xticks([1 101 201 301 401])
% xticklabels({'-1', '-0.5', '0', '0.5', '1'})
% a = get(gca,'XTickLabel'); set(gca,'xticklabel',a,'fontsize',18);
% xlabel('time relative to go cue')
% subplot(1,2,2);
% plot_cov(N_right_sorted',cutout); title('Right')
% xlabel('time relative to go cue')
% sgtitle('psth heatmap - sorted by mean move epoch activity')
% xticks([1 101 201 301 401])
% xticklabels({'-1', '-0.5', '0', '0.5', '1'})
% a = get(gca,'XTickLabel'); set(gca,'XTickLabel',a,'fontsize',18)
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
% % fig_name = fullfile(fig_pth,whos_data,'psth_heatmap');
% % if sav; saveas(gcf,fig_name,'png'); end

% % cov
% cutout = 0.01;
% figure
% subplot(2,2,1)
% plot_cov(Cprep,cutout); title('Prep Covariance')
% subplot(2,2,2)
% plot_cov(Cmove,cutout); title('Move Covariance')
% subplot(2,2,3)
% imagesc(imresize(Cprep,1/6,'bilinear')); title('Prep Downsampled'); colorbar;
% subplot(2,2,4)
% imagesc(imresize(Cmove,1/6,'bilinear')); title('Move Downsampled'); colorbar;
% fig_name = fullfile(fig_pth,whos_data,'cov');
% % if sav; saveas(gcf,fig_name,'png'); end

%% Null space of Q_move (is it equal to Q_prep)

% % null space of Q_move
% [U,S,V] = svd(Q_move);
% nulldim = 

% full_left_null = full_psth{1}*nulldim';
% full_right_null = full_psth{2}*nulldim';
% figure
% plot(full_time,full_left_null(:,1),'r')
% hold on
% plot(full_time,full_right_null(:,1),'b')



%% Helper Functions

function plot_cov(data,cutout)
    num_to_cut = ceil( numel(data) * cutout / 2);
    sorted_data = sort(data(:));
    cmin = sorted_data( num_to_cut );
    cmax = sorted_data( end - num_to_cut + 1);
    imagesc(data, [cmin, cmax]);
    colorbar;
end % plot_cov

function my_animate(gcf,fig_pth,whos_data)
    axis tight
    set(gca,'xticklabel',[]);
    set(gca,'yticklabel',[]);
    set(gca,'zticklabel',[]);
    az = 0;
    el = 90;
    view([az,el]);
    degStep = 1;
    detlaT = 0.1;
    fCount = 71;
    f = getframe(gcf);
    [im,map] = rgb2ind(f.cdata,256,'nodither');
    im(1,1,1,fCount) = 0;
    k = 1;
    % spin 45°
    for i = 0:-degStep:-45
      az = i;
      view([az,el]);
      f = getframe(gcf);
      im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
      k = k + 1;
    end
    % tilt down
    for i = 90:-degStep:15
      el = i;
      view([az,el]);
      f = getframe(gcf);
      im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
      k = k + 1;
    end
    % spin left
    for i = az:-degStep:-90
      az = i;
      view([az,el]);
      f = getframe(gcf);
      im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
      k = k + 1;
    end
    % spin right
    for i = az:degStep:0
      az = i;
      view([az,el]);
      f = getframe(gcf);
      im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
      k = k + 1;
    end
    % tilt up to original
    for i = el:degStep:90
      el = i;
      view([az,el]);
      f = getframe(gcf);
      im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
      k = k + 1;
    end
    fig_name = fullfile(fig_pth,whos_data,'state_space.gif');
    imwrite(im,map,fig_name,'DelayTime',detlaT,'LoopCount',inf)
end % my_animate

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
        plot(time,left(:,i),'r','LineWidth',2)
        hold on
        plot(time,right(:,i),'b','LineWidth',2); xlim(x_lim); ylim(y_lim);
        set(gca,'yticklabel',[]);
        if i~=d; set(gca,'xticklabel',[]); 
        else; set(gca,'fontsize',18); xlabel('time relative to go cue'); end
        hold off
        sgtitle(title_str)
    end
end