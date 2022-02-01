function rez = subspaceid_optimization_moveLabeled(obj, params)

%% PREPROCESS
rez.psth = obj.psth(:,:,params.condToUse);
rez.time = obj.time;
rez.condition = params.condition(params.condToUse);
rez.trialid = params.trialid(params.condToUse);

rez = preprocess_optimization(rez); % softnorma and mean center across conditions

%% GET COVARIANCE MATRICES of EPOCHS
rez = epochCovariance(rez);

%% OPTIMIZATION

% find number of dims for move and prep epochs
rez = epochDimensionality(rez,params.varToExplain);

% main optimization step
rez.alpha = 0; % regularization hyperparam (+ve->discourage sparity, -ve->encourage sparsity)
[Q,~,P,~,~] = orthogonal_subspaces(rez.Cmove,rez.dMove, ... 
                                    rez.Cprep,rez.dPrep,rez.alpha);
                                

rez.Qpotent = Q*P{1};
rez.Qnull = Q*P{2};

rez = getVarianceExplained(rez);

%% PLOTS
% 
% cols = params.cols;
% plotLatents(rez,obj,meta,params, 'Optimization');
% lbl = {'Potent 1', 'Potent 2', 'Null 1'};
% condLbl = {meta.condition{params.conditions}};
% plotStateSpaceGUI(obj.time, obj.psth, rez, cols, 'Optimization', params.dims, lbl, condLbl);
% cond = params.conditions;
% plotSingleTrialsGUI(obj.time,obj.trialpsth,rez,cols,'Optimization', params.dims,lbl,condLbl,meta.trialNum,cond);

clrs = {[0 0 1],[1 0 0],[0.5 0.5 1],[1 0.5 0.5],'k'};

Q = rez.Qpotent;

for i = 1:size(Q,2) % dimension
    f = figure;
    title(['Dim ' num2str(i) '   |   %VE: ' num2str(rez.Qpotent_ve(i))]);
    xlim([params.tmin, params.tmax]);
    hold on
    for j = 1:4%size(obj.psth,3) % condition
        proj = obj.psth(:,:,j) * Q(:,i);
        plot(obj.time, mySmooth(proj,100), 'Color', clrs{j}, ...
            'LineWidth', 2.5);
    end
    hold off
    ax = gca;
    ax.FontSize = 20;
end


end % subspaceid_optimization

%% Helper Functions

function [psthprep, psthmove] = epochPSTH(rez)
% concatenates psth to (ct x n)
psthprep = rez.psth(rez.prepix,:,1);
psthmove = rez.psth(rez.moveix,:,1);
for i = 2:size(rez.psth,3)
    psthprep = [psthprep ; rez.psth(rez.prepix,:,i)];
    psthmove = [psthmove ; rez.psth(rez.moveix,:,i)];
end
end % epochPSTH

function rez = epochCovariance(rez)
% computes covariance of each epoch
[psthprep, psthmove] = epochPSTH(rez);
rez.Cprep = cov(psthprep);
rez.Cmove = cov(psthmove);
end % epochCovariance

function rez = epochDimensionality(rez, varToExplain)
    [psthprep, psthmove] = epochPSTH(rez);
    
    [~,~,explained] = myPCA(psthprep);
    rez.dPrep = numComponentsToExplainVariance(explained, varToExplain);
        
    [~,~,explained] = myPCA(psthmove);
    rez.dMove = numComponentsToExplainVariance(explained, varToExplain);
    
    rez.dMax = max(rez.dMove,rez.dPrep);
end % epochDimensionality


function rez = getVarianceExplained(rez)
prepeigs = sort(eig(rez.Cprep),'descend');
moveeigs = sort(eig(rez.Cmove),'descend');

prepproj = rez.Qnull'*rez.Cprep*rez.Qnull;
moveproj = rez.Qpotent'*rez.Cmove*rez.Qpotent;

rez.Qnull_ve = diag(prepproj) / sum(prepeigs) * 100;
rez.Qpotent_ve = diag(moveproj) / sum(moveeigs) * 100;

end % getVarianceExplained









