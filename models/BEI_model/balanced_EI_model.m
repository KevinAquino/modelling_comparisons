% Function here to simulate the balanaced Excitation-Inhibition model
% 
% This calculates 1 instance and you can feed in previously calculated J values
% i.e. the balancing terms. This is there because this takes too long to do each
% time. 
function bds = balanced_EI_model(C,we,J,total_time,TR)
    if(nargin<3)
        disp('Balancing the model - can take a long time');
        J=Balance_J(we,C);
    end

    N=size(C,2);

    dtt = 1e-3;   % Sampling rate of simulated neuronal activity (seconds)
    dt=0.1;

    load_BEI_parameters;

    % Add an extra 10 seconds to the simulation for equilibration
    Tmaxneuronal=(total_time+10)/dtt;


    
    disp([' At global coupling strength: G=' num2str(we)]);
    neuro_act=zeros(Tmaxneuronal,N);
    sn=0.001*ones(N,1);
    sg=0.001*ones(N,1);
    nn=1;
    for t=0:dt:Tmaxneuronal
        xn=I0*Jexte+w*JN*sn+we*JN*C*sn-J.*sg;
        xg=I0*Jexti+JN*sn-sg;
        rn=phie(xn);
        rg=phii(xg);
        sn=sn+dt*(-sn/taon+(1-sn)*gamma_factor.*rn./1000.)+sqrt(dt)*sigma_factor*randn(N,1);
        sn(sn>1) = 1;
        sn(sn<0) = 0;
        sg=sg+dt*(-sg/taog+rg./1000.)+sqrt(dt)*sigma_factor*randn(N,1);
        sg(sg>1) = 1;
        sg(sg<0) = 0;        
        if abs(mod(t,1))<0.01
            neuro_act(nn,:)=rn';
            nn=nn+1;
        end
    end
    nn=nn-1;

    % %%%% BOLD empirical
    % Friston BALLOON MODEL
    T = nn*dtt; % Total time in seconds

    B = BOLD(T,neuro_act(1:nn,1)'); % B=BOLD activity, bf=Foutrier transform, f=frequency range)
    BOLD_act = zeros(length(B),N);
    BOLD_act(:,1) = B;

    for nnew=2:N
        B = BOLD(T,neuro_act(1:nn,nnew));
        BOLD_act(:,nnew) = B;
    end

    BOLD_sampling=floor(TR/1e-3);
    bds=BOLD_act(1:BOLD_sampling:end,:);

    % Now retrive the total volumes, removing the first 5
    nFrames=floor(total_time/TR);
    bds=bds(5:end,:);
    % keyboard % Change this entirely

    % % Total time for the neural activity
    % time = 0:dtt:nn*dtt;
    % % Conversion to ms
    % time = time*1e3;
    % % Now downsample to 0.1s
    % dsRate=100;
    % % Inputs to the BOLD model
    % z=neuro_act(1:dsRate:end,:);
    % time=time(1:dsRate:end);

    % % Here change it so we use the BOLD_model instead

    % for n=1:size(C,1),
    %     bold_estimate(n,:) = BOLD_model(z(:,n),time,TR);
    % end

    % bds=bold_estimate(:,5:end)';