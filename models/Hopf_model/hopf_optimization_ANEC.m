% This code runs the hopf model and optimizes the weighted connectivity matrix,
% to find the optimal weighted matrix.

% In practise this is done with 80% of the data and tested on 20% (run multiple times?)
% For starters maybe dont have to do this.

function [ts_simulated,Coptim,MSE_run] = hopf_optimization_ANEC(C,G,f_diff,time_series,TR,total_time,GSR_flag)


    % Grabbing all the parameters from the time series
    [Tmax,N,NSUB] = size(time_series);
    Tmax=Tmax*NSUB; % This just ensures a longer time period for better estimation of the phases

    % Calculate the peak frequency - this is done outside the function
    omega = repmat(2*pi*f_diff',1,2); omega(:,1) = -omega(:,1);

    % Setting of the simulation parameters
    dt=0.1*TR/2;
    sig=0.04; % This sigma is used to determine the strength of noise in the Hopf model.
    a=-0.01*ones(N,2);

    % Calculation of the instantenous phase for the time series (this is different to the FCD calculation)
    for nsub=1:NSUB
    	% Calculate instantaneous FC (this is a sub-function below)
    	iFC = calculate_iFC(time_series(:,:,nsub)',TR);	
        FCphasesemp2(nsub,:,:)=squeeze(mean(iFC));    
    end
    FCphasesemp=squeeze(mean(FCphasesemp2));

    % This step here now runs the bulk of the code, we also run this for every weight independently
    % (does this need to be done all at once as a starting point - maybe.)
    %%%%%%%%%%%%
    %% Optimize
    %%
    WE=G;
    Cnew=C;
    % ITER=1:30;
    ITER=[];

    iwe=1;
    for we=WE,        
    	% disp(we);
        Cnew=C;
        % Could possibly use parfor here for parallel processing.
        MSE=100;    
        MSE_run = [];
        wC = we*Cnew;
        for iter=ITER
            disp(['G=',num2str(we),', ITERATION_=',num2str(iter),' MSE=',num2str(MSE)]);
            wC = we*Cnew;
            sumC = repmat(sum(wC,2),1,2); % for sum Cij*xj
            if length(ITER)>1
                % keyboard
            	% Solve the hopf ode/sde:

            	% Solve the HOPF model here.
                % keyboard
            	xs = solve_hopf_ode(omega,a,wC,dt,Tmax,TR,sig);
            	nn = size(xs,2);
                %%%%
                BOLD=xs';
                % Now calculate the instantenous FC using the phase type of analysis
                % keyboard
                iFCsim = calculate_iFC(BOLD,TR);

                FCphases=squeeze(mean(iFCsim));
                
                % Implementation of ANEC
                for i=1:N
                    for j=i+1:N
                        if (C(i,j)>0 || j==N/2+i)
                            Cnew(i,j)=Cnew(i,j)+0.01*(FCphasesemp(i,j)-FCphases(i,j));
                            if Cnew(i,j)<0
                                Cnew(i,j)=0;
                            end
                            Cnew(j,i)=Cnew(i,j);
                        end
                    end
                end
                
                Cnew=Cnew/max(max(Cnew))*0.2;
                
                D = abs(FCphasesemp-FCphases).^2;            
                MSE = sum(D(:))/numel(FCphases);
                MSE_run = [MSE_run,MSE];
                if MSE<0.01
                    break;
                end
            end
        end
        
        Coptim(iwe,:,:)=Cnew;

        % Added part here was not correct.
        wC=we*Cnew;
        
        %%%%%%%%%%%%%%
        %%% Final simulation to save for further use.
        xs = solve_hopf_ode(omega,a,wC,dt,Tmax,TR,sig);    
        
        iwe=iwe+1;
        ts_simulated(:,:,iwe) = xs;
    end

end

% A sub-function to calculate instantneous FC
function iFC = calculate_iFC(time_series,TR)
    % keyboard
	delt = TR;            % sampling interval
	k=2;                  % 2nd order butterworth filter
	fnq=1/(2*delt);
	flp = .04;           % lowpass frequency of filter
	fhi = .07;           % highpass
	Wn=[flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
	[bfilt2,afilt2]=butter(k,Wn);   % construct the filter
	[N,Tmax] = size(time_series);

	signaldata=squeeze(time_series);  %%seed, time, sub    
    Phase_BOLD_data=zeros(N,Tmax);
    % keyboard
    for seed=1:N
        signaldata(seed,:)=signaldata(seed,:)-mean(signaldata(seed,:));
        signal_filt_data =filtfilt(bfilt2,afilt2,signaldata(seed,:));
        Phase_BOLD_data(seed,:) = angle(hilbert(signal_filt_data));
    end
    % Pre-allocation makes this 100 times faster
    iFC=zeros(Tmax,N,N);

    for t=1:Tmax
        for n=1:N
            for p=1:N
                iFC(t,n,p)=cos(Phase_BOLD_data(n,t)-Phase_BOLD_data(p,t));
            end
        end
    end

end