function [states_ordered, params_ordered] = order_states(states,params)
%%%%%%%%%%%%%%%%%%%%
% Orders the states from an HMM-Bayes run based on their parameters.
% D states (no V) come first, ordered from largest to smallest D value
% (smaller D may indicate more complex behavior, e.g. confinement).
% DV states come next, ordered from smallest to largest V magnitude.
%
%%%%%%%%%%%%%%%%%%%%
% Copyright MIT 2015
% Laboratory for Computational Biology & Biophysics
%%%%%%%%%%%%%%%%%%%%


% D states (no V) come first, ordered from largest to smallest D value
% (smaller D may indicate more complex behavior, e.g. confinement)
D_states = find(sum(params.mu_emit,1)==0);
[~,D_idx] = sort(params.sigma_emit(D_states),'descend');
final_order = D_states(D_idx);

% DV states come next, ordered from smallest to largest V magnitude
DV_states = find(sum(params.mu_emit,1)~=0);
V_mag = sum(params.mu_emit(:,DV_states).*params.mu_emit(:,DV_states),1);
[~,DV_idx] = sort(V_mag);
final_order = [final_order DV_states(DV_idx)];


% Re-order the states
if iscell(states)
    states_ordered = cell(size(states));
    for i=1:length(states)
        states_ordered{i} = states{i} + max(states{i});
        for j=1:max(states{i})
            states_ordered{i}(states_ordered{i}==final_order(j)+max(states{i})) = j;
        end
    end
else
    states_ordered = states + max(states);
    for i=1:max(states)
        states_ordered(states_ordered==final_order(i)+max(states)) = i;
    end
end


% Re-order the parameters
params_ordered = params;
params_ordered.p_start = params.p_start(final_order);
params_ordered.p_trans = params.p_trans(final_order,final_order);
params_ordered.mu_emit = params.mu_emit(:,final_order);
params_ordered.sigma_emit = params.sigma_emit(final_order);


end
