clear all;

Data=load('BOLDData_Petra');
indexsub=find(Data.AGES<35);
insub=1;
for nsub=indexsub
    SCdata(insub,:,:)=Data.(sprintf('SUBJECT_%d', nsub)).SC;
    insub=insub+1;
end
C=squeeze(mean(SCdata,1));
C=C/max(max(C))*0.2;

load  empiricalICschizos.mat;

load APARC_parcellation_data.mat;

N=68;
Tmax=148;
Isubdiag = find(tril(ones(N),-1));

% %%%%%%%%%%%%%%
CONTROLS=1:128;
SCHIZOS=129:186;

CASESUB=CONTROLS;
we= ;


NSUB=length(CASESUB);
Projectiondata=Projection1;
clear Projection1;

% %%%%%%%%%%%%%%

TR=2.;  % Repetition Time (seconds)

delt = TR;            % sampling interval
k=2;                  % 2nd order butterworth filter
fnq=1/(2*delt);
flp = .04;           % lowpass frequency of filter
fhi = .07;           % highpass
Wn=[flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
[bfilt2,afilt2]=butter(k,Wn);   % construct the filter


%%%%%%%%%%%%%%%%%%


for nsub=CASESUB
    clear PowSpect;
    signaldata=squeeze(ICA_AROMA_2P_GSR(:,:,nsub));
    TT=Tmax;
    Ts = TT*TR;
    freq = (0:TT/2-1)/Ts;
    [aux minfreq]=min(abs(freq-0.04));
    [aux maxfreq]=min(abs(freq-0.07));
    nfreqs=length(freq);
    
    for seed=1:N
        x=detrend(demean(signaldata(seed,:)));
        ts =zscore(filtfilt(bfilt2,afilt2,x));
        pw = abs(fft(ts));
        PowSpect(:,seed,insub) = pw(1:floor(TT/2)).^2/(TT/TR);
    end
    insub=insub+1;
end

Power_Areas=mean(PowSpect,3);
for seed=1:N
    Power_Areas(:,seed)=gaussfilt(freq,Power_Areas(:,seed)',0.01);
end

[maxpowdata,index]=max(Power_Areas);
f_diff = freq(index);

clear PowSpect  Power_Areas ;

%%%%%%%%%%%%%%%%%%

omega = repmat(2*pi*f_diff',1,2); omega(:,1) = -omega(:,1);

dt=0.1*TR/2;
sig=0.02;
dsig = sqrt(dt)*sig; % to avoid sqrt(dt) at each time step

%%%%%%%%%%%%
%% Optimize
%%

%%%%%%%%%%%%%%

a=zeros(N,2);
ITER=1:100;

nsub2=1;

for nsub=CASESUB
    
    signaldata=squeeze(ICA_AROMA_2P_GSR(:,:,nsub));
    Phase_BOLD_data=zeros(N,Tmax);
    for seed=1:N
        signaldata(seed,:)=signaldata(seed,:)-mean(signaldata(seed,:));
        signal_filt_data =filtfilt(bfilt2,afilt2,signaldata(seed,:));
        Phase_BOLD_data(seed,:) = angle(hilbert(signal_filt_data));
    end
    
    for t=1:Tmax
        for n=1:N
            for p=1:N
                iFC(t,n,p)=cos(Phase_BOLD_data(n,t)-Phase_BOLD_data(p,t));
            end
        end
    end
    FCphasesemp=squeeze(mean(iFC));
    
    Cnew=C;
    Tmax1=Tmax;
    for iter=ITER
        wC = we*Cnew;
        sumC = repmat(sum(wC,2),1,2); % for sum Cij*xj
        xs=zeros(Tmax1,N);
        %number of iterations, 100 willk�hrlich, weil reicht in diesem Fall
        z = 0.1*ones(N,2); % --> x = z(:,1), y = z(:,2)
        nn=0;
        % discard first 3000 time steps
        for t=0:dt:3000
            suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
            zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
            z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(N,2);
        end
        % actual modeling (x=BOLD signal (Interpretation), y some other oscillation)
        for t=0:dt:((Tmax1-1)*TR)
            suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
            zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
            z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(N,2);
            if abs(mod(t,TR))<0.01
                nn=nn+1;
                xs(nn,:)=z(:,1)';
            end
        end
        
        %%%%
        BOLD=xs';
        signal_filt=zeros(N,nn);
        Phase_BOLD=zeros(N,nn);
        
        for seed=1:N
            BOLD(seed,:)=BOLD(seed,:)-mean(BOLD(seed,:));
            signal_filt =filtfilt(bfilt2,afilt2,BOLD(seed,:));
            Phase_BOLD(seed,:) = angle(hilbert(signal_filt));
        end
        
        for t=1:size(BOLD,2)
            for n=1:N
                for p=1:N
                    iFCsim(t,n,p)=cos(Phase_BOLD(n,t)-Phase_BOLD(p,t));
                end
            end
        end
        FCphases=squeeze(mean(iFCsim));
        
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
        if MSE<0.01
            break;
        end   
    end
    Ceff(nsub2,:,:)=Cnew;
    nsub2=nsub2+1;
end



save effconn_control_schizos.mat Ceff;


