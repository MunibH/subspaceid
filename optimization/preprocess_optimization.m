function [obj, meta] = preprocess_optimization(meta, obj)

%% soft-normalize 
meta.lambda = 0.01;
psth = softnorm(obj.psth, meta.lambda);

%% mean-center
psth = meancenter(psth);

%% 
obj.psth = psth;

end % preprocess_optimization

%% 
function snpsth = softnorm(psth, lambda)
% firing rate of a neuron, x of size (time,1), is transformed as:
% x_norm = x / (lambda + max(x) - min(x))
snpsth = psth ./ (lambda + max(psth) - min(psth));
end % softnorm

function mcpsth = meancenter(psth)
% compute mean activity of each neuron across conditions at each time point
% subtract mean from each condition's response

for clu = 1:size(psth,2)
    % find mean at each time point for each condition (time,1)
    mean_across_cond = mean(psth(:,clu,:),3);
    psth(:,clu,:) = psth(:,clu,:) - repmat(mean_across_cond,[1,1,2]);
end
mcpsth = psth;

end % meancenter







