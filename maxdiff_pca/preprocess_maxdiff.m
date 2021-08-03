function [obj, meta] = preprocess_maxdiff(meta, obj)

%% soft-normalize 
meta.lambda = 0.01;
psth = softnorm(obj.psth, meta.lambda);

%% mean-center
psth = meancenter(psth);

%% 
obj.psth = psth;

end % preprocess_maxdiff

%% 
function snpsth = softnorm(psth, lambda)
% firing rate of a neuron, x of size (time,1), is transformed as:
% x_norm = x / (lambda + max(x) - min(x))
snpsth = psth ./ (lambda + max(psth) - min(psth));
end % softnorm

function mcpsth = meancenter(psth)
mcpsth = psth - mean(psth);
end % meancenter