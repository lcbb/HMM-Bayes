function states = hmm_viterbi(obs,p_start,p_trans,mu_emit,sigma_emit)
%%%%%%%%%%%%%%%%%%%%
% Performs the Viterbi algorithm to calculate the most likely sequence of
% hidden states given the HMM parameters and the sequence of observations,
% assuming the emission probability distributions are Gaussian.
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
pointers = zeros(K,T-1);
states = zeros(1,T);

% First time point: F(k) = log(p(y1|k)) + log(pi(k))
for k=1:K
    F(k,1) = -log((2*pi)^(d/2)*sigma_emit(k)^d) - sum((obs(:,1)-mu_emit(:,k)).^2)/(2*sigma_emit(k)^2) + log(p_start(k));
end

% Other time points: F(k,t) = log(p(yt|k)) + max_k'(log(phi(k',k))+F(k',t-1))
for t=2:T
    
    % Recursion
    for k=1:K
        [max_prev pointers(k,t-1)] = max(F(:,t-1)+log(p_trans(:,k)));
        F(k,t) = -log((2*pi)^(d/2)*sigma_emit(k)^d) - sum((obs(:,t)-mu_emit(:,k)).^2)/(2*sigma_emit(k)^2) + max_prev;
    end
    
end

% Final maximization: max_k(F(k,T))
[~, states(T)] = max(F(:,T));

% Reconstruct sequence of states from pointers
for t=T-1:-1:1
    states(t) = pointers(states(t+1),t);
end


end
