function errors = hmm_param_errors(samples_mcmc,Vstates)
%%%%%%%%%%%%%%%%%%%%
% Calculates standard deviations in parameter values from an MCMC
% run (samples_mcmc).
%
%%%%%%%%%%%%%%%%%%%%
% Copyright MIT 2015
% Laboratory for Computational Biology & Biophysics
%%%%%%%%%%%%%%%%%%%%


[d, K] = size(samples_mcmc(1).mu_emit);
Vnonzero = find(Vstates~=0);

errors = struct;


for p = 1:K  % starting probabilities
    
    list = struct2vector(samples_mcmc,'p_start',p);
    errors.p_start(p) = std(list);
  
end


for col = 1:K  % transition probabilities

    for row = 1:K

        list = struct2vector(samples_mcmc,'p_trans',row,col);
        errors.p_trans(row,col) = std(list);
        
    end

end


errors.mu_emit = zeros(d,K);

for col = Vnonzero  % emission means

    for row = 1:d

        list = struct2vector(samples_mcmc,'mu_emit',row,col);
        errors.mu_emit(row,col) = std(list);
        
    end

end


for p = 1:K  % standard deviations
    
    list = struct2vector(samples_mcmc,'sigma_emit',p);
    errors.sigma_emit(p) = std(list);
  
end


end