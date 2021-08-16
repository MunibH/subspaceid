function rez = subspaceid_optimization(meta, obj, params)

%% PREPROCESS
obj.psth = obj.psth(:,:,params.conditions);
rez.time = obj.time;
rez.psth = obj.psth;
rez.trialpsth = obj.trialpsth;
rez.condition = meta.condition(params.conditions);

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
                                
% P1 = [eye(rez.dMove); zeros(rez.dPrep,rez.dMove)];
% P2 = [zeros(rez.dMove, rez.dPrep); eye(rez.dPrep)];

rez.Qpotent = Q*P{1};
rez.Qnull = Q*P{2};

rez.Qpotent_ve = ones(rez.dMove,1);
rez.Qnull_ve   = ones(rez.dPrep,1);

%% PLOTS

cols = {[0,0,1],[1,0,0]};
plotLatents(obj.time, obj.psth, rez, meta, cols, 'Optimization');
% plotStateSpace(obj.time, obj.psth, rez, cols, 'Optimization', params.dims);
lbl = {'Potent 1', 'Potent 2', 'Null 1'};
condLbl = {meta.condition{params.conditions}};
plotStateSpaceGUI(obj.time, obj.psth, rez, cols, 'Optimization', params.dims, lbl, condLbl);
cond = params.conditions;
plotSingleTrialsGUI(obj.time,obj.trialpsth,rez,cols,'Maxdiff', params.dims,lbl,condLbl,meta.trialNum,cond);

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
    
    [~,explained] = myPCA(psthprep);
    rez.dPrep = numComponentsToExplainVariance(explained, varToExplain);
        
    [~,explained] = myPCA(psthmove);
    rez.dMove = numComponentsToExplainVariance(explained, varToExplain);
    
    rez.dMax = max(rez.dMove,rez.dPrep);
end % epochDimensionality


