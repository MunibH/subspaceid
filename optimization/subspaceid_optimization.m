function rez = subspaceid_optimization(meta, obj, params)

%% PREPROCESS
obj.psth = obj.psth(:,:,params.conditions);
rez.time = obj.time;
rez.psth = obj.psth;
rez.trialpsth = obj.trialpsth;
rez.condition = meta.condition(params.conditions);
rez.trials = meta.trialNum(params.conditions);

[obj, meta] = preprocess_optimization(meta, obj);

%% PREP and MOVE EPOCHS
% prepix and moveix corresponds to time idx for each epoch
[rez.prepix, rez.moveix] = getEpochs(obj.time, meta.prepEpoch, meta.moveEpoch);

%% GET COVARIANCE MATRICES of EPOCHS
rez = epochCovariance(obj.psth, rez);

%% OPTIMIZATION

% find number of dims for move and prep epochs
rez = epochDimensionality(rez, obj.psth, rez.prepix, rez.moveix, params.varToExplain);

% main optimization step
rez.alpha = 0; % regularization hyperparam (+ve->discourage sparity, -ve->encourage sparsity)
[Q,~,P,~,~] = orthogonal_subspaces(rez.Cmove,rez.dMove, ... 
                                    rez.Cprep,rez.dPrep,rez.alpha);
                                

rez.Qpotent = Q*P{1};
rez.Qnull = Q*P{2};

rez = getVarianceExplained(rez);

%% PLOTS

cols = {[0,0,1],[1,0,0]};
plotLatents(obj.time, obj.psth, rez, meta, cols, 'Optimization');
lbl = {'Potent 1', 'Potent 2', 'Null 1'};
condLbl = {meta.condition{params.conditions}};
plotStateSpaceGUI(obj.time, obj.psth, rez, cols, 'Optimization', params.dims, lbl, condLbl);
cond = params.conditions;
plotSingleTrialsGUI(obj.time,obj.trialpsth,rez,cols,'Optimization', params.dims,lbl,condLbl,meta.trialNum,cond);

end % subspaceid_optimization

%% Helper Functions

function [psthprep, psthmove] = epochPSTH(psth, prepix, moveix)
% concatenates psth to (ct x n)
psthprep = psth(prepix,:,1);
psthmove = psth(moveix,:,1);
for i = 2:size(psth,3)
    psthprep = [psthprep ; psth(prepix,:,i)];
    psthmove = [psthmove ; psth(moveix,:,i)];
end
end % epochPSTH

function rez = epochCovariance(psth, rez)
% computes covariance of each epoch
[psthprep, psthmove] = epochPSTH(psth, rez.prepix, rez.moveix);
rez.Cprep = cov(psthprep);
rez.Cmove = cov(psthmove);
end % epochCovariance

function rez = epochDimensionality(rez, psth, prepix, moveix, varToExplain)
    [psthprep, psthmove] = epochPSTH(psth, prepix, moveix);
    
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









