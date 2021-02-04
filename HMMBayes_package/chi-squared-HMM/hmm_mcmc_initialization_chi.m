function [best_sample mcmc_params samples logprobs] = hmm_mcmc_initialization_chi(obs,d,K,mu_max,sigma_max,mcmc_params)
%%%%%%%%%%%%%%%%%%%%
% Tests a series of random initializations of MCMC sampling of the
% parameter space of an HMM, to find the best starting parameters
%
% obs - 1xT vector of observations (squared displacement magnitudes)
%       (or a cell of such observations)
% d - number of dimensions in original data
% K - number of states in the HMM
% mu_max - initial mu values will be drawn from [0 mu_max]
%          (the mean parameter for a chi distribution is the magnitude
%           of the vector of means of its gaussian components)
% sigma_max - initial sigma values will be drawn from [0 sigma_max] 
%
% mcmc_params:
%    .Vstates - length K vector of whether each state has zero or nonzero V (e.g. [0 0 1] for K=3 with one nonzero V state)
%    .proposaltype - 'uniform' or 'gaussian' (form of the proposal distribution)
%    .move - 'all' or 'single' or 'block' (number of parameters to move at each step)
%    .delta_p_start - amount by which starting probabilities can be moved in a single step
%    .delta_p_trans - amount by which transition probabilities can be moved in a single step
%    .delta_mus - amount by which emission probability means can be moved in a single step
%    .delta_sigmas - amount by which emission probability standard deviations can be moved in a single step
%    .nTrials - (optional) number of MCMC initialization restarts
%
%%%%%%%%%%%%%%%%%%%%
% Copyright MIT 2015
% Laboratory for Computational Biology & Biophysics
%%%%%%%%%%%%%%%%%%%%


if isfield(mcmc_params,'nTrials')
    nTrials = mcmc_params.nTrials;
elseif K>1
    nTrials = 20;
else
    nTrials = 5;
end


mcmc_params.nIter = 5000;
mcmc_params.delta_p_start = 0.2;
mcmc_params.delta_p_trans = 0.05;

mcmc_params.delta_mus = mu_max/50;
mcmc_params.delta_sigmas = sigma_max/20;

p_start = ones(1,K) / K;
p_trans = ones(K) / K;

%success = zeros(1,nTrials);
best_logprob = -Inf;
best_sample = [];

samples = cell(1,nTrials);
logprobs = cell(1,nTrials);

for i=1:nTrials
    
    mu_emit = rand(1,K) * mu_max;
    mu_emit(mcmc_params.Vstates==0) = 0;  % enforce states with zero V
    
    sigma_emit = rand(1,K) * sigma_max;

    [samples{i} logprobs{i} accept_rate] = hmm_mcmc_chi(obs,d,p_start,p_trans,mu_emit,sigma_emit,mcmc_params);
    
    %figure, plot(logprobs{i})
    
    [max_logprob idx] = max(logprobs{i});

    if max_logprob > best_logprob 
        best_logprob = max_logprob;
        best_sample = samples{i}(idx);
    end
    
    if ~isempty(find(mcmc_params.Vstates~=0,1))
        if accept_rate.mu_emit < 0.3 || accept_rate.mu_emit > 0.5
            mcmc_params.delta_mus = mcmc_params.delta_mus * exp(accept_rate.mu_emit/0.4 - 1.25);
        end
    end
    
    if accept_rate.sigma_emit < 0.3 || accept_rate.sigma_emit > 0.5
        mcmc_params.delta_sigmas = mcmc_params.delta_sigmas * exp(accept_rate.sigma_emit/0.4 - 1.25);
    end
    
    
end


end