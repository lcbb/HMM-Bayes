function hmm_skeleton()
% HMM_SKELETON An example file demonstrating how to create a matrix of
% displacements from a trajectory, initialize HMM parameters, run the HMM
% on that series of displacements, and then plot the result.

    % Load a trajectory

    data = importdata('data/example_track.mat');
    
    track = data.track_fig1;
    cfg = data.cfg_fig1;

    % cfg is a structure with three members:
    % cfg.umperpx -> micrometers per pixel conversion factor
    % cfg.fs -> framerate of the movie
    % cfg.locerror -> localization error
    
    % Set up parameters for the HMM-Bayes algorithm
    
    mcmc_params.parallel = 'on'; % turn off if parallel is not available
    maxK = 3; % maximum number of hidden states to test

    % Create displacements from a [# dimensions x track length] trajectory
    
    steps = track(:,2:end) - track(:,1:end-1);
    
    % Run the algorithm
    % Compile the results into a structure for feeding into the plotting
    % routine and later analysis.
    
    [results.PrM, results.ML_states, results.ML_params, results.full_results, full_fitting, results.logI]...
        = hmm_process_dataset(steps,maxK,mcmc_params);
    
    results.track = track;
    results.steps = steps;
    
    % Save the results of the analysis
    %
    % full_fitting is stored in a separate .mat because it can be quite
    % large and results in long .mat load times. This file contains the
    % samples of model parameters from Markov Chain Monte Carlo.
    
    save('data/analysis_output','cfg','results')
    save('data/analysis_output_samples','cfg','results','full_fitting')

    % Run the plotting algorithm
    
    hmm_results_plot(cfg,results);
end