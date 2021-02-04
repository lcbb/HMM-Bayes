function logI = hmm_integration_gaussian(obs,samples_mcmc,N,Vstates)
%%%%%%%%%%%%%%%%%%%%
% Calculates the numerical integral of the likelihood of the HMM
% observations by Monte Carlo sampling, where the sampling distribution is
% given by a Gaussian approximation for each mu_emit and sigma_emit
% (using the mean and standard deviaiton of the MCMC samples) and uniform
% sampling for the probabilities (from a uniform simplex).
%
% Number of Gaussian parameters:
%   K' * d  : emission probability means (where K' is the number of nonzero V states)
%   K       : emission probability standard deviations
%
%%%%%%%%%%%%%%%%%%%%
% Copyright MIT 2015
% Laboratory for Computational Biology & Biophysics
%%%%%%%%%%%%%%%%%%%%


%tic

[d K] = size(samples_mcmc(1).mu_emit);
Vnonzero = find(Vstates~=0);

% Find the Gaussian approximation for the MCMC sampled distributions of
% mu_emit and sigma_emit

mu_emit_means = [];
mu_emit_stdevs = [];
sigma_emit_means = [];
sigma_emit_stdevs = [];

for col = Vnonzero  % emission means

    for row = 1:d

        list = struct2vector(samples_mcmc,'mu_emit',row,col);
        mu_emit_means(end+1) = mean(list);
        mu_emit_stdevs(end+1) = std(list);
        
    end

end

for p = 1:K  % standard deviations
    
    list = struct2vector(samples_mcmc,'sigma_emit',p);
    sigma_emit_means(end+1) = mean(list);
    sigma_emit_stdevs(end+1) = std(list);
  
end



% Generate N samples from the Gaussian approximations

if ~isempty(Vnonzero)
    samples_mu = randn(N,length(Vnonzero)*d);
    samples_mu = bsxfun(@times,samples_mu,mu_emit_stdevs);
    samples_mu = bsxfun(@plus,samples_mu,mu_emit_means);
end

samples_sigma = randn(N,K);
samples_sigma = bsxfun(@times,samples_sigma,sigma_emit_stdevs);
samples_sigma = bsxfun(@plus,samples_sigma,sigma_emit_means);
samples_sigma = abs(samples_sigma);



% Calculate integral (expected value of the likelihood divided by the
% sampling distribution, calculated over all samples)

logratios = [];
simplex_pdf = (factorial(K-1)/sqrt(K))^(K+1);  % K+1 simplexes are sampled in each iteration below

for i=1:N
    
    p_start = sample_simplex(K);
    
    for j = 1:K
        p_trans(j,:) = sample_simplex(K);
    end
    
    mu_emit = zeros(d,K);
    if ~isempty(Vnonzero)
        mu_emit(:,Vnonzero) = reshape(samples_mu(i,:),d,length(Vnonzero));
    end
    sigma_emit = samples_sigma(i,:);
    
    
    % Get the likelihood of the observations at this sampling point
    if iscell(obs)
        logprob = 0;
        for j=1:length(obs)
            logprob = logprob + hmm_forward(obs{j},p_start,p_trans,mu_emit,sigma_emit);
        end
    else
        logprob = hmm_forward(obs,p_start,p_trans,mu_emit,sigma_emit);
    end
    
    
    % Find the sampling distribution pdf at this sampled point:
    
    if ~isempty(Vnonzero)
        mu_pdf = 1/((2*pi)^(length(Vnonzero)*d/2)*prod(mu_emit_stdevs)) * exp(-sum(((samples_mu(i,:)-mu_emit_means)./mu_emit_stdevs).^2)/2);
    else
        mu_pdf = 1;
    end
    
    for j = 1:K  % sigma pdf is f(x)+f(-x) (because of abs function when sampling sigma above)
        sigma_pdf(j) = 1/(sqrt(2*pi)*sigma_emit_stdevs(j)) * (exp(-(((samples_sigma(i,j)-sigma_emit_means(j))/sigma_emit_stdevs(j))^2)/2) + exp(-(((-samples_sigma(i,j)-sigma_emit_means(j))/sigma_emit_stdevs(j))^2)/2));
    end
    
    sample_pdf = simplex_pdf * mu_pdf * prod(sigma_pdf);
    
    
    % Ratio of likelihood to sampling pdf
    logratios(i) = logprob - log(sample_pdf);
    
end

logratios = logratios(~isnan(logratios));

%logI = log(mean(exp(logratios)));   %this expression can cause underflow
R = max(logratios);
logI = R + log(mean(exp(logratios-R)));   %trick to avoid underflow


% Uniform prior on parameters

prior_scale = 200;  % multiple of the std dev in each parameter

if ~isempty(Vnonzero)
    prior = simplex_pdf * 1/(prod(2*prior_scale*mu_emit_stdevs)*prod(2*prior_scale*sigma_emit_stdevs));
else
    prior = simplex_pdf * 1/prod(2*prior_scale*sigma_emit_stdevs);
end

logI = log(prior) + logI;


%toc

end



