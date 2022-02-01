function rez = preprocess_maxdiff(rez)

%% soft-normalize 
rez.softnorm_lambda = 0.01;
psth = softnorm(rez.psth, rez.softnorm_lambda);

%% mean-center
psth = meancenter(psth);

%% 
rez.psth = psth;

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