function [Joptim p1]=Balance_J_analytic(we,C)

	load_BEI_parameters;
	% keyboard
	% A rescaling here to make it easier.
	taon=taon/1000;
	taog=taog/1000;



	% Solved by working out the transcendtal equation:
	ii=linspace(0,1,10000);	
	functional=Jexti*I0 + JN*gamma_factor*r_0*taon/(1+ gamma_factor*r_0*taon) - taog*phii(ii) - ii;
	[~,ind]=min(abs(functional));
	ii=ii(ind);

	disp(['Solved transcendtal equation numerically, Ii(I)=',num2str(ii)]);

	Joptim = sum(C)'*we*JN/(phii(ii)*taog)*(gamma_factor*r_0*taon/(1+ taon*gamma_factor*r_0)) + 1;
