function [PrM, ML_states, ML_params, full_results, full_fitting] = hmm_process_dataset_chi(steps,Kmax,mcmc_params)
%%%%%%%%%%%%%%%%%%%%
% Performs chi-squared HMM-Bayes analysis of a timeseries of displacements
% (steps) from a particle trajectory.
% If steps is a cell, containing steps from multiple independent
% trajectories, the likelihood is calculated as a product over the
% individual likelihoods.
%
% steps - dxT vector/matrix of observations (for d-dimensional data)
%         (or a cell of such observations)
% Kmax - maximum number of states in the HMM to test
% mcmc_params.parallel - 'on' or 'off' for parallel processing of models
%
%%%%%%%%%%%%%%%%%%%%
% Copyright MIT 2015
% Laboratory for Computational Biology & Biophysics
%%%%%%%%%%%%%%%%%%%%


tStart = tic;

% Ensure that steps has dimensions 1xn where n is the number of tracks.
if iscell(steps)
    steps = reshape(steps,1,numel(steps));
end

% Remove steps with NaN values
if iscell(steps)
    for i=1:length(steps)
        steps{i} = steps{i}(:,~isnan(sum(steps{i},1)));
    end
else
    steps = steps(:,~isnan(sum(steps,1)));
end

% Remove extraneous dimension if present
if iscell(steps)
    for i=1:length(steps)
        if isempty(find(steps{i}(end,:)~=0,1))
            steps{i} = steps{i}(1:end-1,:);
        end
    end
else
    if isempty(find(steps(end,:)>0,1))
        steps = steps(1:end-1,:);
    end
end

% Choose ranges for MCMC parameter initializations
if iscell(steps)
    steps_matrix = cell2mat(steps);
else
    steps_matrix = steps;
end
mu_max_dim = max(abs(steps_matrix),[],2); % maximum mu in each dimension
mu_max = norm(mu_max_dim);
%sigma_max = max(std(steps_matrix,[],2)); % maximum sigma in all dimensions
sigma_max = mu_max;

nModels = sum((1:Kmax)+1); % # models (with and without V states for each K)


Klist = cell(0);
Vstates_list = cell(0);
for K=1:Kmax
    for i=0:K % for each of the K+1 models with K states
        Klist{end+1} = K;
        Vstates_list{end+1} = [zeros(1,K-i) ones(1,i)]; % Vstates vector for each of the models
    end
end

logI = zeros(1,nModels);
full_results = struct('PrM',[],'K',Klist,'Vstates',Vstates_list,'ML_states',[],'ML_params',[]);
full_fitting = struct('mcmc_params_final',[],'logI',cell(1,nModels),'logprobs',[],'samples',[]);

if ~isfield(mcmc_params,'proposaltype')
    mcmc_params.proposaltype = 'gaussian';
end
if ~isfield(mcmc_params,'move')
    mcmc_params.move = 'block';
end


% Convert steps from vectors to squared magnitudes to use as observations
% in the chi-distribution HMM (squared displacement magnitudes)
if iscell(steps)
    for i=1:length(steps)
        obs{i} = sum(steps{i}.^2,1);
    end
    d = size(steps{1},1);
else
    obs = sum(steps.^2,1);
    d = size(steps,1);
end

% Find the marginal likelihood of the data for each number of states
if isfield(mcmc_params,'parallel') && strcmp(mcmc_params.parallel,'on')
    
    if matlabpool('size') == 0, matlabpool; end
    
    parfor m = 1:nModels
        
        K = Klist{m};
        full_fitting(m).mcmc_params_final = mcmc_params;
        full_fitting(m).mcmc_params_final.Vstates = Vstates_list{m};

        % Run initializations with different starting parameters
        [best_sample, full_fitting(m).mcmc_params_final] = hmm_mcmc_initialization_chi(obs,d,K,mu_max,sigma_max,full_fitting(m).mcmc_params_final);

        % Run a longer MCMC trajectory with the best starting parameters
        full_fitting(m).mcmc_params_final.nIter = 100000;
        [full_fitting(m).samples, full_fitting(m).logprobs] = hmm_mcmc_chi(obs,d,best_sample.p_start,best_sample.p_trans,best_sample.mu_emit,best_sample.sigma_emit,full_fitting(m).mcmc_params_final);

        % Perform integration to get marginal likelihood
        nIntegration = 200000;
        logI(m) = hmm_integration_gaussian_chi(obs,d,full_fitting(m).samples(ceil(full_fitting(m).mcmc_params_final.nIter/2):end),nIntegration,full_fitting(m).mcmc_params_final.Vstates);
        full_fitting(m).logI = logI(m);
        
    end
    
    
else
    
    for m = 4 %1:nModels

        K = Klist{m};
        full_fitting(m).mcmc_params_final = mcmc_params;
        full_fitting(m).mcmc_params_final.Vstates = Vstates_list{m};

        % Run initializations with different starting parameters
        [best_sample, full_fitting(m).mcmc_params_final] = hmm_mcmc_initialization_chi(obs,d,K,mu_max,sigma_max,full_fitting(m).mcmc_params_final);

        % Run a longer MCMC trajectory with the best starting parameters
        full_fitting(m).mcmc_params_final.nIter = 100000;
        [full_fitting(m).samples, full_fitting(m).logprobs] = hmm_mcmc_chi(obs,d,best_sample.p_start,best_sample.p_trans,best_sample.mu_emit,best_sample.sigma_emit,full_fitting(m).mcmc_params_final);

        % Perform integration to get marginal likelihood
        nIntegration = 200000;
        logI(m) = hmm_integration_gaussian_chi(obs,d,full_fitting(m).samples(ceil(full_fitting(m).mcmc_params_final.nIter/2):end),nIntegration,full_fitting(m).mcmc_params_final.Vstates);
        full_fitting(m).logI = logI(m);
        
    end
end

PrM = exp(logI-max(logI))/sum(exp(logI-max(logI)));  %subtract max to avoid underflow


for m=4 %1:nModels
    
    full_results(m).PrM = PrM(m);
    
    % Maximum likelihood parameters for each model
    [~, ML_idx] = max(full_fitting(m).logprobs);
    full_results(m).ML_params = full_fitting(m).samples(ML_idx);
    
    % Viterbi algorithm to get most probable state sequence
    if iscell(steps)
        full_results(m).ML_states = cell(size(steps));
        for i=1:length(steps)
            full_results(m).ML_states{i} = hmm_viterbi_chi(obs{i},d,full_results(m).ML_params.p_start,full_results(m).ML_params.p_trans,full_results(m).ML_params.mu_emit,full_results(m).ML_params.sigma_emit);
        end
    else
        full_results(m).ML_states = hmm_viterbi_chi(obs,d,full_results(m).ML_params.p_start,full_results(m).ML_params.p_trans,full_results(m).ML_params.mu_emit,full_results(m).ML_params.sigma_emit);
    end
    
    % Order D states from lowest to highest D and then DV states from lowest to highest V magnitude
    [full_results(m).ML_states, full_results(m).ML_params] = order_states(full_results(m).ML_states, full_results(m).ML_params);
    
end

% Parameters and states of the best model
[~, ML_M] = max(PrM);
ML_states = full_results(ML_M).ML_states;
ML_params = full_results(ML_M).ML_params;


toc(tStart)


end

