% % http://pillowlab.princeton.edu/teaching/statneuro2018/slides/notes03a_SVDandLinSys.pdf
function null_potent_regression_alm_irn()
    clear,clc,close all
    % load data
    data_pth = '/Users/Munib/Documents/Economo-Lab/subspace-id/data/kaufman/';
    whos_data = 'erin';
    data_file = 'EEL6_2020-03-01_processed_alm.mat'; 
    load(fullfile(data_pth,whos_data,data_file));
    
    % reduce dim of irn and alm data
    [~,alm] = pca([data.psth{1}(:,data.alm_idx),data.psth{2}(:,data.alm_idx)]);
    [~,irn] = pca([data.psth{1}(:,data.irn_idx),data.psth{2}(:,data.irn_idx)]);
    alm = alm(:,1:6);
    irn = irn(:,1:6);
    
    % cross validate to find regularization parameter to use
    lambdas = linspace(-1e5,1e5,1000);
%     lambda = cross_validate(irn(data.move_idx,:),alm(data.move_idx,:),lambdas);
    lambda = 10000;
    % compute transformation matrix, W, using ridge regression
    W = ridge_regression(irn(data.move_idx,:),alm(data.move_idx,:),lambda);
    
    % find null and potent spaces of W (reference above)
    [U,S,V] = svd(W);
    W_potent = U(logical(diag(S)),:);
    W_null = U(~sum(S'),:);
    % The number of null and potent dimensions isn't 3 and 3 with this
    % data. Check eigenvalues
    W_potent = U(1:5,:);
    W_null = U(6,:);
    
    plot_projections(data,data.time,W_null,'Null')
    plot_projections(data,data.time,W_potent,'Potent')
    
    % % prep tuning
    numDims = 10;
    gamma = norm((W_null * alm(data.move_idx,:)'),'fro')^2 / ...
            norm((W_potent * alm(data.move_idx,:)'),'fro')^2;
    
    N_null = (W_null * alm(data.prep_idx,:)')';
    N_potent = (W_potent * alm(data.prep_idx,:)')';
    
    tuning_ratio = (1/gamma) * (norm(N_null,'fro')^2 / norm(N_potent,'fro')^2)
    'hi'
end % null_potent_regression


%% Regression and Cross Validation Functions

function W = ridge_regression(Y,X,lambda)
% returns transformation matrix W
% inputs:
% Y - target matrix (T,M)
% X - data matrix   (T,N)
% lambda = regularization parameter

    [~,M] = size(Y);
    [~,N] = size(X);

    Xplus = [X; sqrt(lambda)*eye(N)];
    Yplus = [Y; zeros(N,M)];

    W = Xplus\Yplus;

end % ridge_regression

function best_lambda = cross_validate(Y,X,lambdas)
    n = size(X,1); % num observations (time)
    K = 5; % fraction of data held out for testing
    c = cvpartition(n,'Kfold',K);
    % for each alpha
    %   for each k in k-fold
    %       keep fold k as hold-out data
    %       use remaining folds and current alpha to estimate W
    %       predict hold-out data: M_test,k = X_test,k * W
    %       compute MSE: |M - M_test,k|^2
    %   end for k
    %   average MSE over the k folds: 1/K * sum(MSE)
    % end for alpha
    % choose optimal value: alpha_opt = argmin_p (MSE_p)
    
    mse_i = zeros(1,length(lambdas));
    for i = 1:length(lambdas)
        mse_k = zeros(1,K);
        for k = 1:K
            idxTrain = training(c,k);
            idxTest = test(c,k);
            W_test = ridge_regression(Y(idxTrain,:),X(idxTrain,:),lambdas(i));
            Y_test = X(idxTest,:) * W_test;
            mse_k(k) = immse(Y_test,Y(idxTest,:));
        end
        mse_i(i) = mean(mse_k);
    end
    [~,minIdx] = min(mse_i);
    best_lambda = lambdas(minIdx);
    
    figure
    plot(lambdas,mse_i,'.','MarkerSize',5);   
    xlabel('lambda')
    ylabel('mse')
    title([num2str(K) '-fold' ' cross-validation'])
end

function plot_projections(data,time,W,titlestr)
    [~,alm_l] = pca(data.psth{1}(:,data.alm_idx));
    [~,alm_r] = pca(data.psth{2}(:,data.alm_idx));
    alm_l = alm_l(:,1:6);
    alm_r = alm_r(:,1:6);
    N_l = (W * alm_l')';
    N_r = (W * alm_r')';
    y_lim(1) = min( [min(min(N_l)), min(min(N_r))] );
    y_lim(2) = max( [max(max(N_l)), max(max(N_r))] );
    
    figure
    for i = 1:size(N_l,2)
        subplot(size(N_l,2),1,i)
        plot(time,N_l(:,i),'b','LineWidth',2); ylim(y_lim); xlim([-2,2]); hold on
        plot(time,N_r(:,i),'r','LineWidth',2); hold on
        fill([-1.15,-0.05,-0.05,-1.15],[y_lim(1),y_lim(1),y_lim(2),y_lim(2)],'g','FaceAlpha',0.05,'EdgeColor','none'); % fill prep epoch
        hold on
        fill([0.05,1.15,1.15,0.05],[y_lim(1),y_lim(1),y_lim(2),y_lim(2)],'r','FaceAlpha',0.05,'EdgeColor','none'); % fill move epoch
        set(gca,'yticklabel',[]);
        if i~=size(N_l,2); set(gca,'xticklabel',[]); 
        else; set(gca,'fontsize',18); xlabel('time relative to go cue'); end
        hold off;
        sgtitle(titlestr)
    end
end % plot_projections


