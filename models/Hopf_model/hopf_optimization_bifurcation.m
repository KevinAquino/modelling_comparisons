% Rewriting entirely the code for the Hopf Model
% Removing all model analyses from here and only calculating the model itself.


% Needs the following subfunctions
% 1. Calculation of peak frequency (can be done easily enough)
% 2. Calculation of power terms
% 3. Calculation of Hopf ODE

% Fed into..
%  Greedy search optimization

function [ts_simulated_all,bifpar,all_a,sc] = hopf_optimization_bifurcation(C,G,f_peak,time_series,TR,total_time,GSR_flag)

if(nargin<7)
    % Do GSR by default!
    GSR_flag=1;
end

rng('shuffle');
nNodes = length(C);

max_hopf_iterations = 200;
% Calculating the power spectrum integral
empirical_power_ratio = calculatePowerRatio(time_series,total_time/TR,TR);

vmax=max(empirical_power_ratio); %consider computing this later where needed
vmin=min(empirical_power_ratio);%consider computing this later where needed

% Previously calculated frequency peak
f_diff = f_peak;
omega = repmat(2*pi*f_diff,1,2); %angular velocity
omega(:,1) = -omega(:,1);

a = repmat(-0.05*ones(nNodes,1),1,2);

% Simulation parameters
dt = 0.1;

minm=100;
sig=0.04;
abort_critera = 0.1;
delta_update = 0.1;
% bestIter = 1; %JVS for tracking purposes

for iter = 1:max_hopf_iterations,
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Adding extra sim time
    xs = solve_hopf_ode(omega,a,G*C,dt,total_time/TR*10,TR,sig);            
    if(GSR_flag)
        xs = RegressNoiseSignal(xs.',mean(xs.')).';
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    sim_power_ratio = calculatePowerRatio(xs,size(xs,1),TR);
    
    vsmin=min(sim_power_ratio);%empirical_power_ratio y sim_power_ratio max min points not equal
    vsmax=max(sim_power_ratio);
    bb=(vmax-vmin)/(vsmax-vsmin);
    aa=vmin-bb*vsmin;  %% adaptation of local bif parameters
    sim_power_ratio=aa+bb*sim_power_ratio;
    minm1=max(abs(empirical_power_ratio-sim_power_ratio)./empirical_power_ratio);    
    
    if minm1<minm
        minm=minm1;
        a1=a;
        bestIter = iter; %JVS: tracking
        best_sim_power_ratio = sim_power_ratio; %JVS: tracking
    end
    
    %--------------------------------------------------------------
    %FEEDBACK
    % %--------------------------------------------------------------
    if(~mod(iter,50))
        % fprintf(1, 'iter: %03d, minm1: %5.3f\n', iter, minm1);
        disp(['iter:',num2str(iter),' , minm1: ',num2str(minm1),', corr: ',num2str(corr(a(:,1),empirical_power_ratio'))])
    end
    sc(iter)=minm1;
    all_a(:,:,iter) = a;
    if(iter>999)
        disp('! Max number of Iterations, did not reach stopping critera');
        disp('Choosing the best');
        smoothed_obj_function = movmean(sc,20);
        [~,index_min] = min(smoothed_obj_function);
        a = all_a(:,:,index_min);
        % Force a finish
        minm=0.05;
    end
    %--------------------------------------------------------------
    
    %CRITERION REACHED?
    if minm<abort_critera,%default is 0.1
        break;
    end
    
    %UPDATE a VALUES FOR NEXT ITER
    if ~any(isnan(sim_power_ratio))
        a(:,1)=a(:,1)+delta_update*(1-sim_power_ratio./empirical_power_ratio)';
        a(:,2)=a(:,2)+delta_update*(1-sim_power_ratio./empirical_power_ratio)';
    else
        %JVS: this should not happen according to Gustavo, but it
        %can when a values run crazy
        warning('There are NaNs in the power spectra. Probably a-values too strong.');
        a = a1;
    end
end
a=a1;
%%%%%
%Linearize for all Gs
% aLinear=aLin(a(:,1)',Cfg.opt_a.gref,wG);

% Below is done for reference
a=a1;%use those avalues that have been found to be optimal
bifpar=a(:,1)';

ts_simulated_all = xs;
  
end


