clear all;

%FIND PATH
path0=pwd;

%%ADD REQUIRED SUBFLODERS
path2=[path0 '/toolbox'];
addpath(path2)

load APARC_parcellation_data.mat;
N=68;
Isubdiag = find(tril(ones(N),-1));

% %%%%%%%%%%%%%%
CONTROLS=1:128; %110;
n_subjects_Controls=length(CONTROLS);
SCHIZOS=129:186;
n_subjects_Schizos=length(SCHIZOS);

nsub1=1;
for nsub=CONTROLS
    signaldata=(squeeze(ICA_AROMA_2P_GSR(:,:,nsub)));
    tc_aal{nsub1,1}=signaldata;
    nsub1=nsub1+1;
end
nsub1=1;
for nsub=SCHIZOS
    signaldata=(squeeze(ICA_AROMA_2P_GSR(:,:,nsub)));
    tc_aal{nsub1,2}=signaldata;
    nsub1=nsub1+1;
end

n_Task=2;
[N_areas, Tmax]=size(tc_aal{1,1});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1 - Compute the ICs
disp('Processing the ICs from BOLD data')

Tmaxtotal=Tmax*(n_subjects_Schizos+n_subjects_Controls);
Tmaxcontrol=Tmax*n_subjects_Controls;
Tmaxschizos=Tmax*n_subjects_Schizos;

% Parameters of the data
TR=2.;  % Repetition Time (seconds)

% Preallocate variables to save FC patterns and associated information
phaselock_all=zeros(Tmaxtotal,length(Isubdiag)); % All leading eigenvectors
Time_all=zeros(2, Tmaxtotal); % vector with subject nr and task at each t
t_all=0; % Index of time (starts at 0 and will be updated until n_Sub*Tmax)

% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.04;                    % lowpass frequency of filter (Hz)
fhi = 0.07;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter
clear fnq flp fhi Wn k

n_Subjects(1)=n_subjects_Controls;
n_Subjects(2)=n_subjects_Schizos;

for task=1:n_Task
    for s=1:n_Subjects(task)
        % Get the BOLD signals from this subject in this task
        [N_areas, Tmax]=size(tc_aal{s,task});
        % Get the BOLD signals from this subject in this task
        BOLD = tc_aal{s,task};
        BOLD=BOLD(:,1:Tmax);
        Phase_BOLD=zeros(N_areas,Tmax);
        
        for seed=1:N_areas
            BOLD(seed,:)=BOLD(seed,:)-mean(BOLD(seed,:));
            signal_filt =filtfilt(bfilt,afilt,BOLD(seed,:));
            Phase_BOLD(seed,:) = angle(hilbert(signal_filt));
        end
       
        for t=1:Tmax
            iPH=zeros(N_areas,N_areas);
            %Calculate the Instantaneous FC (BOLD Phase Synchrony)
            for n=1:N_areas
                for m=1:N_areas
                    iPH(n,m)=cos(Phase_BOLD(n,t)-Phase_BOLD(m,t));
                end
            end
            t_all=t_all+1; % Update time
            phaselock_all(t_all,:)=iPH(Isubdiag);
            Time_all(:,t_all)=[s task]; % Information that at t_all, V1 corresponds to subject s in a given task
        end
    end
end

ICcomp = assembly_patterns_nozscore(phaselock_all');
NumAssemblies=size(ICcomp,2);

%% 2 - Projections

disp('Projections')

for t=1:Tmaxtotal
    for ass=1:NumAssemblies
        Projection(t,ass)=(dot(ICcomp(:,ass),phaselock_all(t,:)'))^2;%
    end
end
    
disp('%%%%% IC_prob SUCCESSFULLY COMPLETED %%%%%%%')
Projection1=[];
Projection2=[];

edges=[];
delta=max(max(Projection))/60;
x=0;
for i=1:60
    edges=[edges x];
    x=x+delta;
end
            
for task=1:n_Task   % 1, Baselineline, 2, baby face task
    for s=1:n_Subjects(task)
        % Select the time points representing this subject and task
        T=((Time_all(1,:)==s)+(Time_all(2,:)==task))>1;
        
        for ass=1:NumAssemblies      
            if task==1
                [Prob1sub(ass,s,:) xx]=histcounts(Projection(T,ass),edges,'Normalization','probability');
            end
            if task==2
                [Prob2sub(ass,s,:) xx]=histcounts(Projection(T,ass),edges,'Normalization','probability');
            end
        end
        
        if task==1
            Projection1=vertcat(Projection1,Projection(T,:));
        end
        
        if task==2
            Projection2=vertcat(Projection2,Projection(T,:));
        end        
    end
end


for ass=1:NumAssemblies
    [aux pval(ass)]=kstest2(Projection1(:,ass),Projection2(:,ass));
end

cindex=FDR_benjHoch(pval,0.05);

for ass=1:NumAssemblies
    for task=1:n_Task   % 1, Baselineline, 2, baby face task
        for s=1:n_Subjects(task)
            if task==1
                i=find(Prob1sub(ass,s,:));
                Entropy1_sub(ass,s)=-sum(Prob1sub(ass,s,i).*log(Prob1sub(ass,s,i)));
            end
            if task==2
                i=find(Prob2sub(ass,s,:));
                Entropy2_sub(ass,s)=-sum(Prob2sub(ass,s,i).*log(Prob2sub(ass,s,i)));
            end
        end
    end
end

%[pvalentro h]=ranksum(sum(Entropy1_sub),sum(Entropy2_sub));

 a=sum(Entropy1_sub);
 b=sum(Entropy2_sub);
 stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],5000,0.05,'ttest');
 pvalentro=min(stats.pvals);

save empiricalICschizos.mat ICcomp NumAssemblies Entropy1_sub Entropy2_sub Projection1 Projection2;

figure
bar([mean(sum(Entropy1_sub)) mean(sum(Entropy2_sub))],'EdgeColor','w','FaceColor',[.5 .5 .5])
hold on
% Error bar containing the standard error of the mean
errorbar([mean(sum(Entropy1_sub)) mean(sum(Entropy2_sub))],[std(sum(Entropy1_sub))/sqrt(numel(sum(Entropy1_sub))) std(sum(Entropy2_sub))/sqrt(numel(sum(Entropy2_sub)))],'LineStyle','none','Color','k')
set(gca,'XTickLabel',{'Health', 'Schizos'})
if pvalentro<0.05
    plot(1.5,max([mean(sum(Entropy1_sub))+std(sum(Entropy1_sub))/sqrt(numel(sum(Entropy1_sub))) mean(sum(Entropy2_sub))+std(sum(Entropy2_sub))/sqrt(numel(sum(Entropy2_sub)))])+0.1,'*k')
end
ylabel('Entropy');


