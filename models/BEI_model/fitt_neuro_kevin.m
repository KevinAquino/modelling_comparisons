clear all;
load Human_66.mat C FC_emp;
C = (C.' + C)/2;
C=C/max(max(C))*0.2;
N=66;
Isubdiag = find(tril(ones(N),-1));



dtt = 1e-3;   % Sampling rate of simulated neuronal activity (seconds)
dt=0.1;

taon=100;
taog=10;
gamma=0.641;
sigma=0.001;
JN=0.15;
I0=0.382;
Jexte=1.;
Jexti=0.7;
w=1.4;

Tmaxneuronal=100000;  %%%(Tmax+10)*2000;
% WE=0:0.05:5;
% WE=2:1:4;

ii=1;
% for we=WE
J=Balance_J(we,C);
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

FCsimul=corrcoef(bds);






