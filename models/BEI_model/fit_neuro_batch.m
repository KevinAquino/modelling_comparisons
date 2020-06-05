% clear all;
% load Human_66.mat C FC_emp;
% load PetraData.mat C FC_emp;

C = (C.' + C)/2;
C=C/max(max(C))*0.2;
% N=68;
N=size(C,2);
Isubdiag = find(tril(ones(N),-1));
G = linspace(0.1,5,20);
for ind=1:length(G),
    we = G(ind);
    % we = 3.1;
    % we = str2double(getenv('we'));
    %%%%%%%%%%%%%%
    %
    % TT=Tmax;
    % Ts = TT*2;
    % freq = (0:TT/2-1)/Ts;
    % [aux minfreq]=min(abs(freq-0.04));
    % [aux maxfreq]=min(abs(freq-0.07));
    % nfreqs=length(freq);
    %
    %
    % delt = 2;            % sampling interval
    % k=2;                  % 2nd order butterworth filter
    % fnq=1/(2*delt);       % Nyquist frequency
    % flp = .04;           % lowpass frequency of filter
    % fhi = .07;           % highpass
    % Wn=[flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
    % [bfilt2,afilt2]=butter(k,Wn);   % construct the filter
    %
    % %%%%%%%%%%%%%%
    %
    % kk3=1;
    % kk4=1;
    % for nsub=1:NSUB
    %     signaldata=eval(sprintf('TC%d', nsub));
    %     FCe(nsub,:,:)=corrcoef(signaldata');
    %     for seed=1:N
    %         x=demean(detrend(signaldata(seed,:)));
    %         timeseriedata = filtfilt(bfilt2,afilt2,x);    % zero phase filter the data
    %         Xanalytic = hilbert(demean(timeseriedata));
    %         Phases(seed,:) = angle(Xanalytic);
    %     end
    %
    %     T=1:Tmax;
    %
    %     for t=T
    %         kudata=sum(complex(cos(Phases(:,t)),sin(Phases(:,t))))/N;
    %         syncdata(t)=abs(kudata);
    %         %%%%
    %         for i=1:N
    %             for j=1:i-1
    %                 patt(i,j)=cos(adif(Phases(i,t),Phases(j,t)));
    %                 phijdata(i,j,kk4)=patt(i,j);
    %             end
    %         end
    %         pattern(t,:)=patt(Isubdiag);
    %         kk4=kk4+1;
    %     end
    %
    %     npattmax=size(pattern,1);
    %     for t=1:npattmax-2
    %         p1=mean(pattern(t:t+2,:));
    %         for t2=t+1:npattmax-2
    %             p2=mean(pattern(t2:t2+2,:));
    %             phfcddata(kk3)=dot(p1,p2)/norm(p1)/norm(p2);
    %             matrixcdc(t,t2)=dot(p1,p2)/norm(p1)/norm(p2);
    %             kk3=kk3+1;
    %         end
    %     end
    %
    %     metastabilitydata2(nsub)=mean(syncdata);
    % end
    % FCemp=squeeze(mean(FCe,1));
    % metastabilitydata=mean(metastabilitydata2);


    %%%%%%%%

    dtt   = 1e-3;   % Sampling rate of simulated neuronal activity (seconds)
    dt=0.1;

    taon=100;
    taog=10;
    gamma=0.641;
    sigma=0.01;
    JN=0.15;
    I0=0.382;
    Jexte=1.;
    Jexti=0.7;
    w=1.4;

    Tmaxneuronal=100000*3;  %%%(Tmax+10)*2000;
    % WE=0:0.05:5;
    % WE=2:1:4;

    % ii=1;
    % for we=WE
    J=Balance_J(we,C);
    % J=Balance_J_analytic(we,C); % Straight analytic approximation not enough. 
    disp([' At global coupling strength: G=' num2str(we)]);
    %     kk3=1;
    %     kk4=1;
    neuro_act=zeros(Tmaxneuronal,N);
    sn=0.001*ones(N,1);
    sg=0.001*ones(N,1);
    nn=1;
    for t=0:dt:Tmaxneuronal
        xn=I0*Jexte+w*JN*sn+we*JN*C*sn-J.*sg;
        xg=I0*Jexti+JN*sn-sg;
        rn=phie(xn);
        rg=phii(xg);
        sn=sn+dt*(-sn/taon+(1-sn)*gamma.*rn./1000.)+sqrt(dt)*sigma*randn(N,1);
        sn(sn>1) = 1;
        sn(sn<0) = 0;
        sg=sg+dt*(-sg/taog+rg./1000.)+sqrt(dt)*sigma*randn(N,1);
        sg(sg>1) = 1;
        sg(sg<0) = 0;
        j=j+1;
        if abs(mod(t,1))<0.01
            neuro_act(nn,:)=rn';
            nn=nn+1;
        end
    end
    nn=nn-1;

    %%%% BOLD empirical
    % Friston BALLOON MODEL
    T = nn*dtt; % Total time in seconds

    B = BOLD(T,neuro_act(1:nn,1)'); % B=BOLD activity, bf=Foutrier transform, f=frequency range)
    BOLD_act = zeros(length(B),N);
    BOLD_act(:,1) = B;

    for nnew=2:N
        B = BOLD(T,neuro_act(1:nn,nnew));
        BOLD_act(:,nnew) = B;
    end

    bds=BOLD_act(2000:2000:end,:);

    % FCsimul=corrcoef(bds);
    ts_simulated_all(:,:,ind) = bds;    
    figure;imagesc(zscore(bds)')
    drawnow;
end
%%%%%%%%%%%%

%         for seed=1:N
%             ts=detrend(demean(bds(1:Tmax,seed)'));
%             tss = filtfilt(bfilt2,afilt2,ts);
%             Xanalytic = hilbert(demean(tss));
%             Phases(seed,:) = angle(Xanalytic);
%         end
%
%         T=1:Tmax;
%
%         for t=T
%             kudata=sum(complex(cos(Phases(:,t)),sin(Phases(:,t))))/N;
%             sync(t)=abs(kudata);
%             %%%%%
%             for i=1:N
%                 for j=1:i-1
%                     patt(i,j)=cos(adif(Phases(i,t),Phases(j,t)));
%                     phij(i,j,kk4)=patt(i,j);
%                 end
%             end
%             pattern(t,:)=patt(Isubdiag);
%             kk4=kk4+1;
%         end
%
%         npattmax=size(pattern,1);
%         for t=1:npattmax-2
%             p1=mean(pattern(t:t+2,:));
%             for t2=t+1:npattmax-2
%                 p2=mean(pattern(t2:t2+2,:));
%                 phfcd(kk3)=dot(p1,p2)/norm(p1)/norm(p2);
%                 kk3=kk3+1;
%             end
%         end
%
%         metastability(ii)=abs(mean(sync)-metastabilitydata);


% cc=corrcoef(FC_emp(Isubdiag),FCsimul(Isubdiag));
% fitting=cc(2);
%     [H(ii),P(ii),ksdist(ii)]=kstest2(phfcddata,phfcd);
%     ksdist

% fileName = ['G_' num2str(we) '.mat'];
% save(fileName,'fitting','FCsimul','we');
%%%%%%%%%%%%

%     ii=ii+1;
% end

% figure
% plot(WE,fitting,'k');






