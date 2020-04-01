% Structural connectivity test:
function [xs] = run_hopf_no_ts(C,G,f_diff,TR)
	% Normalization of C
	C=C/max(C(:))*0.2;
	% Find the number of nodes
	N = size(C,1);
	wC = C*G;
	sig=0.04;
	dt=0.1;
	% TR=2;
	Tmax=1200;


	% Definition of node frequencies & conversion to radians	
	omega = repmat(2*pi*f_diff,1,2); 
	omega(:,1) = -omega(:,1);

	% Defintion of homgeneous a parameter
	a=-0.01*ones(N,2);	

	% Solve the Hopf model, model parameters
	xs = solve_hopf_ode(omega,a,wC,dt,Tmax,TR,sig);
