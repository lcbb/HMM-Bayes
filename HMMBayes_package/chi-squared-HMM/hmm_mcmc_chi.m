function [samples logprobs accept_rate] = hmm_mcmc_chi(obs,d,guess_p_start,guess_p_trans,guess_mu_emit,guess_sigma_emit,mcmc_params)
%%%%%%%%%%%%%%%%%%%%
% Performs Markov Chain Monte Carlo sampling of the parameter space of an
% HMM
%
% obs - 1xT vector of observations (squared displacement magnitudes)
%       (or a cell of such observations)
% d - number of dimensions in original data
% guess_p_start - 1xK vector of starting probabilities for K states
% guess_p_trans - KxK matrix of transition probabilities for K states
% guess_mu_emit - 1xK vector of emission probability means
%               (the mean parameter for a chi distribution is the magnitude
%               of the vector of means of its gaussian components)
% guess_sigma_emit - 1xK vector of emission probability standard deviations
%
% mcmc_params:
%    .nIter - number of MCMC iterations
%    .Vstates - length K vector of whether each state has zero or nonzero V (e.g. [0 0 1] for K=3 with one nonzero V state)
%    .proposaltype - 'uniform' or 'gaussian' (form of the proposal distribution)
%    .move - 'all' or 'single' or 'block' (number of parameters to move at each step)
%    .delta_p_start - amount by which starting probabilities can be moved in a single step
%    .delta_p_trans - amount by which transition probabilities can be moved in a single step
%    .delta_mus - amount by which emission probability means can be moved in a single step
%    .delta_sigmas - amount by which emission probability standard deviations can be moved in a single step
%
%%%%%%%%%%%%%%%%%%%%
% Copyright MIT 2015
% Laboratory for Computational Biology & Biophysics
%%%%%%%%%%%%%%%%%%%%


% Ensure that means are zero for states with zero V
guess_mu_emit(:,mcmc_params.Vstates==0) = 0;


% Start storing all sampled parameters
samples = repmat(struct('p_start',[],'p_trans',[],'mu_emit',[],'sigma_emit',[]),1,mcmc_params.nIter);
samples(1).p_start = guess_p_start;
samples(1).p_trans = guess_p_trans;
samples(1).mu_emit = guess_mu_emit;
samples(1).sigma_emit = guess_sigma_emit;

% Start storing all log likelihoods
logprobs = zeros(1,mcmc_params.nIter);
if iscell(obs)
    for i=1:length(obs)
        logprobs(1) = logprobs(1) + hmm_forward_chi(obs{i},d,samples(1).p_start,samples(1).p_trans,samples(1).mu_emit,samples(1).sigma_emit);
    end
else
    logprobs(1) = hmm_forward_chi(obs,d,samples(1).p_start,samples(1).p_trans,samples(1).mu_emit,samples(1).sigma_emit);
end

% Keep track of acceptance rate
if strcmp(mcmc_params.move,'single')
    accept.p_start = [];
    accept.p_trans = [];
    accept.mu_emit = [];
    accept.sigma_emit = [];
elseif strcmp(mcmc_params.move,'block')
    accept.probs = [];
    accept.mu_emit = [];
    accept.sigma_emit = [];
elseif strcmp(mcmc_params.move,'all')
    accept = ones(mcmc_params.nIter,1);
end


% Set proposal function to uniform or gaussian
if strcmp(mcmc_params.proposaltype,'uniform')
    f = @(x)rand(x)*2-1;
elseif strcmp(mcmc_params.proposaltype,'gaussian')
    f = @randn;
end

%tic

% Loop
for i=2:mcmc_params.nIter
    
    % Step 1: Propose a move
    if strcmp(mcmc_params.move,'single')
        [samples(i) reject selection] = update_single_chi(samples(i-1),mcmc_params,f);
        accept.(selection)(end+1) = 1;
    elseif strcmp(mcmc_params.move,'block')
        [samples(i) reject selection] = update_block_chi(samples(i-1),mcmc_params,f);
        accept.(selection)(end+1) = 1;
    elseif strcmp(mcmc_params.move,'all')
        [samples(i) reject] = update_all_chi(samples(i-1),mcmc_params,f);
    end
  
    
    % Step 2: Accept or reject move using the Metropolis method
    if ~reject
        
%         tic
        
        if iscell(obs)
            for j=1:length(obs)
                logprobs(i) = logprobs(i) + hmm_forward_chi(obs{j},d,samples(i).p_start,samples(i).p_trans,samples(i).mu_emit,samples(i).sigma_emit);
            end
        else
            logprobs(i) = hmm_forward_chi(obs,d,samples(i).p_start,samples(i).p_trans,samples(i).mu_emit,samples(i).sigma_emit);
        end
        
        if isnan(logprobs(i))
            reject = 1;
        elseif logprobs(i) < logprobs(i-1) && log(rand) > logprobs(i) - logprobs(i-1)  % Reject
            reject = 1;
        end
        
%         toc
%         disp(samples(i))
%         disp(logprobs(i))
        
    end

    
    if reject
        samples(i) = samples(i-1);
        logprobs(i) = logprobs(i-1);
        if strcmp(mcmc_params.move,'single') || strcmp(mcmc_params.move,'block')
            accept.(selection)(end) = 0;
        elseif strcmp(mcmc_params.move,'all')
            accept(i) = 0;
        end
    end
    
end

%toc

if strcmp(mcmc_params.move,'single')
    accept_rate.p_start = sum(accept.p_start)/length(accept.p_start);
    accept_rate.p_trans = sum(accept.p_trans)/length(accept.p_trans);
    accept_rate.mu_emit = sum(accept.mu_emit)/length(accept.mu_emit);
    accept_rate.sigma_emit = sum(accept.sigma_emit)/length(accept.sigma_emit);
elseif strcmp(mcmc_params.move,'block')
    accept_rate.probs = sum(accept.probs)/length(accept.probs);
    accept_rate.mu_emit = sum(accept.mu_emit)/length(accept.mu_emit);
    accept_rate.sigma_emit = sum(accept.sigma_emit)/length(accept.sigma_emit);
elseif strcmp(mcmc_params.move,'all')
    accept_rate = sum(accept)/length(accept);
end


end



function [new reject selection] = update_single_chi(old,mcmc_params,f)

reject = 0;

K = length(old.p_start);
Vnonzero = find(mcmc_params.Vstates~=0);

new = old;

% Number of parameters:
%   K       : starting probabilities (K-length vector that sums to 1)
%   K^2     : transition probabilities (KxK matrix where each row sums to 1)
%   K'      : emission probability means (where K' is the number of nonzero V states)
%   K       : emission probability standard deviations
nParams = K + K^2 + length(Vnonzero) + K;

% Randomly choose a parameter to move
p = randi(nParams,1);


% Propose a move for the selected parameter

if p <= K  % move a starting probability
    
    new.p_start(p) = old.p_start(p) + f(1)*mcmc_params.delta_p_start;
    
    % Reject moves outside 0 and 1
    if new.p_start(p)<0 || new.p_start(p)>1
        reject = 1;
    end
    
    % Re-normalize to 1 
    new.p_start = new.p_start/sum(new.p_start);
    
    selection = 'p_start';
    
    
elseif p <= K + K^2  % move a transition probability
    
    p = p - K;
    row = ceil(p/K);
    col = mod(p-1,K)+1;
    
    new.p_trans(row,col) = old.p_trans(row,col) + f(1)*mcmc_params.delta_p_trans;

    % Reject moves outside 0 and 1
    if new.p_trans(row,col)<0 || new.p_trans(row,col)>1
        reject = 1;
    end
    
    % Re-normalize to 1 
    new.p_trans(row,:) = new.p_trans(row,:)/sum(new.p_trans(row,:));
    
    selection = 'p_trans';
    
    
elseif p <= K + K^2 + length(Vnonzero)  % move an emission mean
    
    p = p - (K + K^2);
    
    new.mu_emit(Vnonzero(p)) = old.mu_emit(Vnonzero(p)) + f(1)*mcmc_params.delta_mus;
    
    % Reject moves below 0 (since mean parameter for a chi distribution is a magnitude)
    if new.mu_emit(Vnonzero(p))<0
        reject = 1;
    end
    
    selection = 'mu_emit';
    
    
else  % move an emission standard deviation
    
    p = p - (K + K^2 + length(Vnonzero));
    
    new.sigma_emit(p) = old.sigma_emit(p) + f(1)*mcmc_params.delta_sigmas;
    
    % Reject moves below 0
    if new.sigma_emit(p)<0
        reject = 1;
    end
    
    selection = 'sigma_emit';
    
end


end




function [new reject selection] = update_block_chi(old,mcmc_params,f)

reject = 0;

K = length(old.p_start);
Vnonzero = find(mcmc_params.Vstates~=0);

new = old;

% Number of parameters:
%   K       : starting probabilities (K-length vector that sums to 1)
%   K^2     : transition probabilities (KxK matrix where each row sums to 1)
%   K'      : emission probability means (where K' is the number of nonzero V states)
%   K       : emission probability standard deviations

% Number of blocks = 3 (probabilities, means, standard deviations)

% Randomly choose a block to move
if ~isempty(Vnonzero)
    b = randi(3,1);  % 3 blocks
else
    b = randi(2,1);  % only 2 blocks (no means to move)
end


% Propose a move for each parameter in the selected block

if b == 1  % move all the probabilities
    
    for p = 1:K  % move all the starting probabilities

        new.p_start(p) = old.p_start(p) + f(1)*mcmc_params.delta_p_start;

        % Reject moves outside 0 and 1
        if new.p_start(p)<0 || new.p_start(p)>1
            reject = 1;
        end

    end

    % Re-normalize to 1 
    new.p_start = new.p_start/sum(new.p_start);


    for row = 1:K  % move all the transition probabilities

        for col = 1:K

            new.p_trans(row,col) = old.p_trans(row,col) + f(1)*mcmc_params.delta_p_trans;

            % Reject moves outside 0 and 1
            if new.p_trans(row,col)<0 || new.p_trans(row,col)>1
                reject = 1;
            end

        end

        % Re-normalize to 1 
        new.p_trans(row,:) = new.p_trans(row,:)/sum(new.p_trans(row,:));

    end
    
    selection = 'probs';
        
    
elseif b == 2  % move all the emission standard deviations
    
    for p = 1:K  

        new.sigma_emit(p) = old.sigma_emit(p) + f(1)*mcmc_params.delta_sigmas;

        % Reject moves below 0
        if new.sigma_emit(p)<0
            reject = 1;
        end

    end

    selection = 'sigma_emit';
    
    
elseif b == 3  % move all the emission means

    for p = 1:length(Vnonzero)  

        new.mu_emit(Vnonzero(p)) = old.mu_emit(Vnonzero(p)) + f(1)*mcmc_params.delta_mus;

        % Reject moves below 0 (since mean parameter for a chi distribution is a magnitude)
        if new.mu_emit(Vnonzero(p))<0
            reject = 1;
        end
        
    end
    
    selection = 'mu_emit';

    
end


end






function [new reject] = update_all_chi(old,mcmc_params,f)

reject = 0;

K = length(old.p_start);
Vnonzero = find(mcmc_params.Vstates~=0);

new = old;

% Number of parameters:
%   K       : starting probabilities (K-length vector that sums to 1)
%   K^2     : transition probabilities (KxK matrix where each row sums to 1)
%   K'      : emission probability means (where K' is the number of nonzero V states)
%   K       : emission probability standard deviations


% Propose a move for each parameter

for p = 1:K  % move all the starting probabilities
    
    new.p_start(p) = old.p_start(p) + f(1)*mcmc_params.delta_p_start;

    % Reject moves outside 0 and 1
    if new.p_start(p)<0 || new.p_start(p)>1
        reject = 1;
    end
    
end

% Re-normalize to 1 
new.p_start = new.p_start/sum(new.p_start);
    

for row = 1:K  % move all the transition probabilities
    
    for col = 1:K
        
        new.p_trans(row,col) = old.p_trans(row,col) + f(1)*mcmc_params.delta_p_trans;
        
        % Reject moves outside 0 and 1
        if new.p_trans(row,col)<0 || new.p_trans(row,col)>1
            reject = 1;
        end
    
    end
    
    % Re-normalize to 1 
    new.p_trans(row,:) = new.p_trans(row,:)/sum(new.p_trans(row,:));
    
end
    
    
for p = 1:length(Vnonzero)  % move all the emission means
    
    new.mu_emit(Vnonzero(p)) = old.mu_emit(Vnonzero(p)) + f(1)*mcmc_params.delta_mus;
    
    % Reject moves below 0 (since mean parameter for a chi distribution is a magnitude)
    if new.mu_emit(Vnonzero(p))<0
        reject = 1;
    end
    
end
    
    
for p = 1:K  % move all the emission standard deviations
    
    new.sigma_emit(p) = old.sigma_emit(p) + f(1)*mcmc_params.delta_sigmas;
    
    % Reject moves below 0
    if new.sigma_emit(p)<0
        reject = 1;
    end
    
end


end



