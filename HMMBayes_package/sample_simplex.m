function p = sample_simplex(n)
%%%%%%%%%%%%%%%%%%%%
% Uniformly samples from an n-dimensional simplex with parameters p_1...p_n
% where: 0<p_i<1 and sum(p_i)=1 (e.g. probabilities)
%
% Algorithm:
%     Set x_0 = 0 and x_n=1.
%     Generate n-1 uniform random draws x_i from the open interval (0,1).
%     Sort into ascending order the n+1 points x_0, ..., x_n.
%     The n coordinates p_1, ..., p_n of the final points on the unit
%     simplex are given by p_i=x_i-x_(i-1). 
% 
% Normalization factor for this uniform sampling pdf is equal to:
%     sqrt(n)/(n-1)!
% ..so to get pdf at the sampled point:
%     prob = factorial(n-1)/sqrt(n);
%
%%%%%%%%%%%%%%%%%%%%
% Copyright MIT 2015
% Laboratory for Computational Biology & Biophysics
%%%%%%%%%%%%%%%%%%%%


x = sort([0; rand(n-1,1); 1]);

% p = circshift(x,-1) - x;
% p = p(1:end-1)';

p = zeros(1,n);
for i=1:n
    p(i) = x(i+1) - x(i);
end


end