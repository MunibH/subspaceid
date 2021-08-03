function psths_for_grant()
% find a few cells from each cell type to show in grant
clear,clc,close all

%% load data
data_pth = '/Users/Munib/Documents/Economo-Lab/subspace-id/data/elsayed/';
whos_data = 'mike';
% data_file = 'alldat.mat'; % all data (includes tagged cells but unidentified)
data_file = 'alldat_tagged_processed_int.mat'; % all data including manually merged tagged cells
load(fullfile(data_pth,whos_data,data_file));

psth = data.psth_tavg;
time = data.full_time;
nontag_idx = data.pop_idx;
ptlow_idx = data.ptlow_idx;
ptup_idx = data.ptup_idx;
clearvars -except psth time nontag_idx ptlow_idx ptup_idx

%% filter to make less squiggly

Fs = 30000;
N = size(psth{1},1);
Nfir = 70;
Fst = 75;
firf = designfilt('lowpassfir','FilterOrder',Nfir, ...
    'CutoffFrequency',Fst,'SampleRate',Fs);
psthf = cellfun(@(x) filter(firf,x),psth,'UniformOutput',false);

% correct time lag
delay = mean(grpdelay(firf));

tt = time(1:end-delay);
psthn = cellfun(@(x) x(1:end-delay,:),psth,'UniformOutput',false);
psthff = psthf;
psthff{1}(1:delay,:) = [];
psthff{2}(1:delay,:) = [];

% plot(tt,psthn{1}(:,1),tt,psthff{1}(:,1)) % check if any time lag


%% plot cells

which_idx = ptlow_idx;
for i = which_idx(1):1:which_idx(end)
    y_lim = [0, max(max(psthff{1}(:,i),psthff{2}(:,i)))];
    plot(tt,psthff{1}(:,i),'r','LineWidth',2); ylim(y_lim); xlim([-2.5,1.5]);
    hold on
    plot(tt,psthff{2}(:,i),'b','LineWidth',2); 
    fill([-1.15,-0.05,-0.05,-1.15],[y_lim(1),y_lim(1),y_lim(2),y_lim(2)],'r','FaceAlpha',0.05,'EdgeColor','none'); % fill prep epoch
    hold on 
    fill([0.05,1.15,1.15,0.05],[y_lim(1),y_lim(1),y_lim(2),y_lim(2)],'g','FaceAlpha',0.05,'EdgeColor','none'); % fill move epoch
    hold off
    title(['Neuron #',num2str(i)])
    pause
end


end % psths_for_grant
