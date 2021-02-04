function hmm_demo_samples()
%HMM_SAMPLES_PLOT Demo of using the full_fitting output variable to plot
%the samples of the fig1 DV model. 

    data = importdata('data/track_fig1_example_samples.mat');
    full_fitting = data.full_fitting;

    % Extract structure variables into vectors

    sigma_samples = struct2vector(full_fitting(2).samples,'sigma_emit',1); 
    mu_x_samples = struct2vector(full_fitting(2).samples,'mu_emit',1);
    mu_y_samples = struct2vector(full_fitting(2).samples,'mu_emit',2);
    % The index 2 above corresponds to the DV model. See Table S1 in the
    % manuscript for the list of fit models in order.
    
    figure(1451)
    clf
    
    subplot(3,1,1)
    hist(sigma_samples,40)
    xlabel('\sigma (A.U.)')
    ylabel('Counts')
    title('Samples of \sigma')
    
    subplot(3,1,2)
    hist(mu_x_samples,40)
    xlabel('\mu_x (A.U.)')
    ylabel('Counts')
    title('Samples of \mu_x')
    
    subplot(3,1,3)
    hist(mu_y_samples,40)
    xlabel('\mu_y (A.U.)')
    ylabel('Counts')
    title('Samples of \mu_y')
    
    suptitle('Markov Chain Monte Carlo Samples of D-V Model Parameters')
    
    set(findall(gcf,'type','text'),'fontSize',12,'fontWeight','bold')
    set(findall(gcf,'type','axes'),'fontSize',14,'fontWeight','bold')
    
    set(gcf,'Position',[0 0 700 900]);
    set(gcf,'color','w');

end

