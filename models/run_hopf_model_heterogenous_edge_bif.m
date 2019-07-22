function [ts_simulated_all,grandFCcorr,FCD,all_anecs,all_MSE] = run_hopf_model_heterogenous_edge_bif(sc_matrix,time_series,G,folder)

% Setting up the parameters
C = sc_matrix;
% N=68;
% Tmax=148;
[N,Tmax] = size(squeeze(time_series(:,:,1)));

Isubdiag = find(tril(ones(N),-1));
NSUB=size(time_series,3);


TR=2.;  % Repetition Time (seconds)

%%%%%%%%%%%%%%
nsub2=1;
phfcddata = [];
for nsub=1:NSUB
	% Calculate instantaneous FC (this is a sub-function below)
	iFC = calculate_iFC(time_series(:,:,nsub),TR);	
    FCphasesemp2(nsub2,:,:)=squeeze(mean(iFC));
    nsub2=nsub2+1;
    FCemp2(nsub2,:,:)=corrcoef(time_series(:,:,nsub).');

    % Add all the phase for fcd together and use it in a universal script. 
    phfcddata = [phfcddata,phase_fcd(time_series(:,:,nsub),TR)];
end

FCphasesemp=squeeze(mean(FCphasesemp2));
FCemp=squeeze(mean(FCemp2));

% Calculate the peak frequency, now this is done in a universal script that is used multiple times.
f_diff = calculate_peak_frequency(time_series,TR);
omega = repmat(2*pi*f_diff',1,2); omega(:,1) = -omega(:,1);

dt=0.1*TR/2;
sig=0.02; % This sigma is used to determine the strength of noise in the Hopf model.

%%%%%%%%%%%%
%% Optimize
%%
WE=G;
a=-0.01*ones(N,2);
Cnew=C;
ITER=1:30;

Tmax=Tmax*NSUB;

iwe=1;
for we=WE,
    figure;
    subplot(1,2,1);
    imagesc(C);
    title('DATA SC');
    axis image;
    disp(['G=',num2str(we)]);
	% disp(we);
    Cnew=C;
    % Could possibly use parfor here for parallel processing.
    MSE=100;    
    MSE_run = [];
    for iter=ITER
        disp(['ITERATION_=',num2str(iter),' MSE=',num2str(MSE)]);
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
    %%% Final simulation.   
    xs = solve_hopf_ode(omega,a,wC,dt,Tmax,TR,sig);    
    
    
    

    % Now save the BOLD time series, as well as the correlation between simulation and experimental signals
    BOLD=xs';      
    ts_simulated_all(:,:,iwe) = BOLD;
	FC_simul=corrcoef(BOLD');
    cc=corrcoef(atanh(FCemp(Isubdiag)),atanh(FC_simul(Isubdiag)));
    grandFCcorr(iwe)=cc(2);

    % Calculate the FCD based on the phases and compare this to the simulation using a KS test implemented in matlab   
    % Note for consistency this is now calculated in a universal function.
    % keyboard
	phfcd = phase_fcd(xs,TR);
    [H,P,FCD(iwe)]=kstest2(phfcd,phfcddata);       

    all_anecs{iwe}=Cnew;
    all_MSE{iwe} = MSE_run;

    subplot(1,2,2);
    imagesc(Cnew);
    title(['ANEC, G =',num2str(we)]);
    axis image;
    drawnow;
    iwe=iwe+1;

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