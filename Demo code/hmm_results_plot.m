function hmm_results_plot(varargin)
% HMM_RESULTS Master function for plotting the results of the HMM analysis.
% Call this function with the output of HMM-Bayes stored in the result
% input variable. cfg stores the parameters of your movie. See
% hmm_skeleton.m for examples of how to create these structures from your
% movie info and the output of the HMM.
%
% INPUTS:
% cfg,result - config and result structures as demonstrated in
% hmm_skeleton.m
%
% OR
%
% mat_file - .mat file you would like to load the results from

    if ischar(varargin{1})
        data = importdata(varargin{1});
        
        cfg = data.cfg;
        result = data.results;
    else
        cfg = varargin{1};
        result = varargin{2};
    end

    % Extract variables from structures

    PrM = result.PrM;
    ML_params = result.ML_params; 

    track = result.track;  
    ML_states = result.ML_states; 
    
    fs = cfg.fs; 
    umperpx = cfg.umperpx; 
    
	% Remove excess NaNs flanking the trajectory
    
    trackidcs = find(~isnan(sum(track,1)));
                
    start_idx = trackidcs(1);
    end_idx = trackidcs(end);
    
    track = track(:,start_idx:end_idx);
    
    steps = result.steps; 

    % Generate a vector used to color the trajectory and the state
    % sequence. This vector contains addresses which serve as indices for
    % the colormap. All that is being done here is accounting for NaNs in
    % the middle of the trajectory - gaps during the tracking process.
    % ML_states contains a vector of numbers ex. 1,2,2,2,1,3,3,2,... which
    % indicate which state the HMM is in at a particular displacement. We
    % use this state sequence for coloring.
    
    color_idcs = ML_states;
    
    for i = 1:size(track,2)-1
        if isnan(sum(sum(track(:,i:i+1))))
            color_idcs = [color_idcs(1:i-1) 0 color_idcs(i:end)];
        end
    end
    
    % Convert track and steps to um
    
    track = track .* umperpx;
    steps = steps .* umperpx;
    
    %%% Draw everything
    
    figure(7141)

    set(gcf,'Position',[0 0 1300 800]);
    clf
    set(gcf, 'color', 'white');
    
    % Generate a color map for coloring the trajectory, state sequence, etc
    % to map pieces back to the states they originated from.
    
    cmap = distinguishable_colors(max(ML_states));
    color_idcs = color_idcs + 1; % Shift everything so indices start at 1 instead of 0.
    
    % Go in sequence drawing the subplots. See individual functions (lower
    % in this file) for descriptions of their function.
    
    subplot(2,2,1); draw_track(track,color_idcs,[0 0 0; cmap]);
    subplot(2,2,2); draw_mlstates(color_idcs,fs,[0 0 0; cmap]);
    subplot(2,4,5); draw_modelprobs(PrM);
    subplot(2,4,6); draw_text(ML_states, ML_params, fs, umperpx, cmap, cfg.locerror)
    subplot(2,4,7); draw_pdf(steps,ML_states,ML_params,umperpx,cmap);
    subplot(2,4,8); draw_stephist(steps,ML_states,cmap);
    
    suptitle('HMM Results')
    
    set(findall(gcf,'type','text'),'fontSize',12,'fontWeight','bold')
    set(findall(gcf,'type','axes'),'fontSize',12,'fontWeight','bold')
end

function draw_track(track,ML_states,cmap)
    %%% Draw the trajectory

    hold all

    % Each point is connected by a colored line whose color is defined by
    % the colormap in the above section. Each color comes from the
    % ML_states vector.
    
    cline(track(1,:),track(2,:),zeros(1,numel(track(1,:))),ML_states,cmap,2);
    
    xmin = min(track(1,:));
    xmax = max(track(1,:));
    ymin = min(track(2,:));
    ymax = max(track(2,:));
    
    xlabel('x (\mum)')
    ylabel('y (\mum)')
    title('Trajectory (ML state sequence)')
    
    axis equal

%     xlim([xmin - 0.5, xmax + 0.5])
%     ylim([ymin - 0.5, ymax + 0.5])

end

function draw_mlstates(ML_states,fs,cmap)
    %%% Draw the ML state sequence

    % One horizontal line is drawn for each timepoint in the Viterbi sequence,
    % which is colored according to the state it falls into.
    
    for i = 1:numel(ML_states)
        cline([i-1 i]/fs,[1 1],[0 0],1,cmap(ML_states(i),:),12);
    end
    
    xlabel('Time (sec)')
    ylabel('State')
    
    set(gca,'YTick',1:max(ML_states))

    ylim([0.5 1.5])
    set(gca, 'TickDir', 'out')
    title('State sequence')
end

function draw_modelprobs(PrM)
    %%% Draw the bar graph for comparison of model probabilities.
    
    bar(PrM)
    
    ylabel('Model probability')
    labels = {'D','DV','D, D','D, DV','DV, DV','D, D, D','D, D, DV','D, DV, DV','DV, DV, DV'};
    set(gca,'XTickLabel',labels)
    rotateXLabels( gca(), 45 )
end

function draw_text(ML_states, ML_params, fs, umperpx, cmap, locerror)
    %%% Calculate average lifetimes of each state. This means the average
    %%% time that a particle spends in each state before switching away.
    
    % How many states do we have?
    uniqstates = unique(ML_states);
    
    % Create a bin for each state that we will be adding times into. These
    % times are the time spent in a particular state at a particular time.
    
    taulist = cell(1,numel(uniqstates));
    
    % Find transition points where we move from one state to another.
    
    transpoints = find(ML_states(1:end-1) - ML_states(2:end))+1;
    
    % Now start adding the lengths of each individual state segment to the
    % bins
    
    curidx = 1;
    
    for i = 1:numel(transpoints)
        state = ML_states(transpoints(i)-1);
        taulist{state} = [taulist{state}, transpoints(i)-curidx];
        
        curidx = transpoints(i);
    end
    
    % Add the length of the last state segment
    
    state = ML_states(end);
    taulist{state} = [taulist{state}, numel(ML_states)-curidx+1];
    
    taubins = zeros(1,numel(uniqstates));
    
    for i = 1:numel(uniqstates)
        taubins(i) = mean(taulist{i});
    end
    
    % Convert from vector length to time
    
    taubins = taubins / fs;
    
    %%% Summary statistics panel - plot the text. The max number of states
    %%% supported here is 3.
    
    if numel(uniqstates) == 1 % One state is present
        hmmlabel={'ML (estimated) parameters:',....
            sprintf('D: %.5f (\\mum^2/s)',sqrt(-2*locerror^2 + ML_params.sigma_emit^2)^2*fs/2*umperpx^2),...
            sprintf('V: [%.5f;%.5f] (\\mum/s)',fs*umperpx*ML_params.mu_emit(1),fs*umperpx*ML_params.mu_emit(2)),...
        };
        text(0,.5,hmmlabel)
    elseif numel(uniqstates) == 2 % Two states are present
        text(0,0.8,'Average lifetimes:');
        text(0,0.7,sprintf('\\tau_1: %.5f (sec)',taubins(1)),'Color',cmap(1,:));
        text(0,0.6,sprintf('\\tau_2: %.5f (sec)',taubins(2)),'Color',cmap(2,:));
        text(0,0.5,'ML (estimated) parameters:');
        text(0,0.3,{sprintf('D_1: %.5f (\\mum^2/s)',sqrt(-2*locerror^2 + ML_params.sigma_emit(1)^2)^2*fs/2*umperpx^2),...
            sprintf('V_1: [%.5f;%.5f] (\\mum/s)',fs*umperpx*ML_params.mu_emit(1,1),fs*umperpx*ML_params.mu_emit(2,1))},'Color',cmap(1,:));
        text(0,0.1,{sprintf('D_2: %.5f (\\mum^2/s)',sqrt(-2*locerror^2 + ML_params.sigma_emit(2)^2)^2*fs/2*umperpx^2),...
            sprintf('V_2: [%.5f;%.5f] (\\mum/s)',fs*umperpx*ML_params.mu_emit(1,2),fs*umperpx*ML_params.mu_emit(2,2))},'Color',cmap(2,:));
    else % Three states are present
        text(0,0.95,'Average lifetimes:');
        text(0,0.85,sprintf('\\tau_1: %.5f (sec)',taubins(1)),'Color',cmap(1,:));
        text(0,0.75,sprintf('\\tau_2: %.5f (sec)',taubins(2)),'Color',cmap(2,:));
        text(0,0.65,sprintf('\\tau_3: %.5f (sec)',taubins(3)),'Color',cmap(3,:));
        text(0,0.55,'ML (estimated) parameters:');
        text(0,0.40,{sprintf('D_1: %.5f (\\mum^2/s)',sqrt(-2*locerror^2 + ML_params.sigma_emit(1)^2)^2*fs/2*umperpx^2),...
            sprintf('V_1: [%.5f;%.5f] (\\mum/s)',fs*umperpx*ML_params.mu_emit(1,1),fs*umperpx*ML_params.mu_emit(2,1))},'Color',cmap(1,:));
        text(0,0.2,{sprintf('D_2: %.5f (\\mum^2/s)',sqrt(-2*locerror^2 + ML_params.sigma_emit(2)^2)^2*fs/2*umperpx^2),...
            sprintf('V_2: [%.5f;%.5f] (\\mum/s)',fs*umperpx*ML_params.mu_emit(1,2),fs*umperpx*ML_params.mu_emit(2,2))},'Color',cmap(2,:));
        text(0,0.0,{sprintf('D_3: %.5f (\\mum^2/s)',sqrt(-2*locerror^2 + ML_params.sigma_emit(3)^2)^2*fs/2*umperpx^2),...
            sprintf('V_3: [%.5f;%.5f] (\\mum/s)',fs*umperpx*ML_params.mu_emit(1,3),fs*umperpx*ML_params.mu_emit(2,3))},'Color',cmap(3,:));
    end
    
    axis off
end

function draw_pdf(steps,ML_states,ML_params,umperpx,cmap)
    %%% Draw the scatterplot of displacements with PDF contours
    
    %%% Generate a scatterplot of the displacements along the trajectory
    
    % Remove NaNs
    steps = steps(:,~isnan(sum(steps)));
    
    hold all
    
    % Plot each set of displacements in the color of the state they
    % correspond to.
    
    for i = unique(ML_states)
        tempsteps = steps(:,i==ML_states);
        plot(tempsteps(1,:),tempsteps(2,:),'.','Color',cmap(i,:))
    end
    
    freezeColors();
    
    % Overlay PDF contour plot on top
    
    mu = ML_params.mu_emit * umperpx;
    sigma = ML_params.sigma_emit * umperpx;

    draw_std(mu,sigma,cmap);
    
    % draw lines showing where 0 is
    
    line(xlim, [0 0],'LineStyle','--','Color','k')
    line([0 0], ylim,'LineStyle','--','Color','k')

    axis equal
    axis tight
    
    xlabel('\Deltax (\mum)')
    ylabel('\Deltay (\mum)')
    title('Trajectory displacements')
    
    freezeColors();
end

function draw_stephist(steps,ML_states,cmap)
    %%% Plot RMS histogram of displacements

    % Calculate length of displacements
    steps_conv = sqrt(sum((steps).^2));
    
    % Remove NaNs
    steps_conv(isnan(steps_conv)) = [];
    
    % Create a histogram of these distances
    
    binranges = linspace(0,max(steps_conv),10);
    [bins,idcs] = histc(steps_conv,binranges);
    
    % Create sub bins based on which state each displacement belongs to
    
    subbins = zeros(numel(bins),numel(unique(ML_states)));
    
    for i = 1:numel(idcs)
        subbins(idcs(i),ML_states(i)) = subbins(idcs(i),ML_states(i)) + 1;
    end
    
    % Plot the histogram
    
    bar(binranges,subbins,1)
    
    hold all

    set(gca, 'TickDir', 'out')
    
    title('Step size distribution')
    xlabel('Step size (\mum)')
    ylabel('Number of steps')
    cur = xlim;
    xlim([0, cur(2)])

    colormap(cmap)
   
end


