function [best_sample, best_mcmc_params, samples, logprobs] = hmm_mcmc_initialization(obs,K,mu_range,sigma_max,mcmc_params)
%%%%%%%%%%%%%%%%%%%%
% Tests a series of random initializations of MCMC sampling of the
% parameter space of an HMM, to find the best starting parameters
%
% obs - dxT vector/matrix of observations (for d-dimensional data)
%       (or a cell of such observations)
% K - number of states in the HMM
% mu_range - initial mu values will be drawn from [-mu_range mu_range]
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
elseif K==1
    nTrials = 5;
elseif K>2 && strfind(which('hmm_forward'),'mex')
    nTrials = 500;
else
    nTrials = 100;
end


d = size(mu_range,1);

mcmc_params.nIter = 5000;
mcmc_params.delta_mus = max(mu_range(:,2)-mu_range(:,1))/50;
mcmc_params.delta_sigmas = sigma_max/50;

p_start = ones(1,K) / K;
p_trans = ones(K) / K;

best_logprob = -Inf;
best_sample = [];
best_mcmc_params = mcmc_params;

samples = cell(1,nTrials);
logprobs = cell(1,nTrials);
mcmc_params_all = cell(1,nTrials+1);
mcmc_params_all{1} = mcmc_params;

for i=1:nTrials
    
    % Generate initial guesses for mu and sigma parameters
    mu_emit = rand(d,K) .* repmat(mu_range(:,2)-mu_range(:,1),1,K) + repmat(mu_range(:,1),1,K);
    mu_emit(:,mcmc_params.Vstates==0) = 0;  % enforce states with zero V
    sigma_emit = rand(1,K) * sigma_max;

    % Run MCMC trial
    [samples{i}, logprobs{i}, accept_rate] = hmm_mcmc(obs,p_start,p_trans,mu_emit,sigma_emit,mcmc_params_all{i});
    
    % Update delta parameters if necessary (but don't allow them to grow larger than the range/4)
    mcmc_params_all{i+1} = mcmc_params_all{i};
    if ~isempty(find(mcmc_params.Vstates~=0,1))
        if accept_rate.mu_emit < 0.3 || accept_rate.mu_emit > 0.5
            mcmc_params_all{i+1}.delta_mus = min(mcmc_params_all{i}.delta_mus*exp(accept_rate.mu_emit/0.4-1.25),max(mu_range(:,2)-mu_range(:,1))/4);
        end
    end
    if accept_rate.sigma_emit < 0.3 || accept_rate.sigma_emit > 0.5
        mcmc_params_all{i+1}.delta_sigmas = min(mcmc_params_all{i}.delta_sigmas*exp(accept_rate.sigma_emit/0.4-1.25),sigma_max/4);
    end
    
    % Store maximum likelihood sample
    [max_logprob, idx] = max(logprobs{i});
    if max_logprob > best_logprob 
        best_logprob = max_logprob;
        best_sample = samples{i}(idx);
        best_mcmc_params = mcmc_params_all{i+1};
    end
    
end


end
