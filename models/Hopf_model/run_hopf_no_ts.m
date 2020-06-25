% Structural connectivity test:
function [xs] = run_hopf_no_ts(C,G,f_diff,TR,total_time,a)	
	disp(['Running the HOPF model for G=',num2str(G)]);
	% Normalization of C
	C=C/max(C(:))*0.2;
	% Find the number of nodes
	N = size(C,1);
	wC = C*G;
	sig=0.04;
	dt=0.1;
	
	% here just adjusting how many volumes one needs for the simulation
	Tmax=ceil(total_time/TR);


	% Definition of node frequencies & conversion to radians	
	omega = repmat(2*pi*f_diff,1,2); 
	omega(:,1) = -omega(:,1);

	% Defintion of homgeneous a parameter
	if(nargin<6)		
		a=-0.01*ones(N,2);	
	end	
	% Solve the Hopf model, model parameters
	xs = solve_hopf_ode(omega,a,wC,dt,Tmax,TR,sig);
