% This function is to precalculate the balanced weights for the DECO 2014 model
% This is done because its easier to do this then recalculating on the fly


% Firstly choose a range of global weights, save them then store them.
% This is done because it is needed for every G.

% Need to put a check within the code to see if this has been pre-calculated.

G = linspace(0,4.5,10);

for g_ind=1:length(G),
	disp(['Precalculation of J_i for G=',num2str(G(g_ind))])
	J_precalculated(g_ind,:) = Balance_J(G(g_ind),C);	
end

save('modes/BEI_model/precalculated_Ji.mat','J_precalculated','G');