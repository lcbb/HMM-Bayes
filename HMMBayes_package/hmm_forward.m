function log_obs_prob = hmm_forward(obs,p_start,p_trans,mu_emit,sigma_emit)
%%%%%%%%%%%%%%%%%%%%
% Performs the forward algorithm to calculate the log probability of the
% sequence of observations given the HMM parameters, assuming the emission
% probability distributions are Gaussian.
%
% obs - dxT vector/matrix of observations (for d-dimensional data)
% p_start - 1xK vector of starting probabilities for K states
% p_trans - KxK matrix of transition probabilities for K states
% mu_emit - dxK vector/matrix of emission probability means
% sigma_emit - 1xK vector of emission probability standard deviations
%
%%%%%%%%%%%%%%%%%%%%
% Copyright MIT 2015
% Laboratory for Computational Biology & Biophysics
%%%%%%%%%%%%%%%%%%%%


T = size(obs,2);
K = length(p_start);
d = size(obs,1);

F = zeros(K,T);
scaling = zeros(1,T);

% First time point: F(k) = p(y1|k)*pi(k)
for k=1:K
    %F(k,1) = mvnpdf(obs(:,1)',mu_emit(:,k)',eye(d)*sigma_emit(k)^2) * p_start(k);
    F(k,1) = 1/((2*pi)^(d/2)*sigma_emit(k)^d) * exp(-sum((obs(:,1)-mu_emit(:,k)).^2)/(2*sigma_emit(k)^2)) * p_start(k);
end
scaling(1) = sum(F(:,1));
F(:,1) = F(:,1)/scaling(1);

% Other time points: F(k,t) = p(yt|k)*sum_k'(phi(k',k)F(k',t-1))
for t=2:T
    
    % Recursion
    for k=1:K
        %F(k,t) = mvnpdf(obs(:,t)',mu_emit(:,k)',eye(d)*sigma_emit(k)^2) * sum(F(:,t-1).*p_trans(:,k));
        F(k,t) = 1/((2*pi)^(d/2)*sigma_emit(k)^d) * exp(-sum((obs(:,t)-mu_emit(:,k)).^2)/(2*sigma_emit(k)^2)) * sum(F(:,t-1).*p_trans(:,k));
    end
    
    % Normalization to 1 (to prevent underflow)
    scaling(t) = sum(F(:,t));
    F(:,t) = F(:,t)/scaling(t);
    
end

% Final probability: sum_k(F(k,T))   %obs_prob = sum(F(:,T));
% Equal to the product of all the scaling factors
log_obs_prob = sum(log(scaling));


end
