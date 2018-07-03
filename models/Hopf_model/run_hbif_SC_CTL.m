%Script to launch the hbif function. Tmax is the number of timepoints, wG is
%the G vector, ldata is the number of subjects and inside the Cfg object
%there are is nIters, which is the number of iterations for initial
%optimization, the update strength of the optimization, the threshold
%abortion criteria (abortCrit) that indicates the maximum error allowed
%between vsig and vsigs, so abortCrit == 0.1 means that if within the
%optimization the error between vsig and vsigs is =< 0.1 the optimization
%stops and uses that last value. If this level is not reached, it uses the
%one with the lowest error. Lastly, gref is the reference G from which
%the initial optimizaiton starts and from which the linearization of
%bifurcation parameters is assumed to start. This might be changed but
%values above 3 are recommended.
%Victor Saenger, 2016

%load('HealthyGroup.mat');%Load data

Tmax = 148;%Check whats your tmax.
wG = 0:0.1:6;%this is de G vector i.e, 0:0.1:8

ldata = 1;%Number of subjets.


Cfg.simulID = 'MortenGroup_v3wide';
Cfg.opt_a.nIters = 200;%Number of iterations for initial optimization
Cfg.opt_a.updateStrength = 0.1;%Update strength for a optimization
Cfg.opt_a.abortCrit = 0.1; %maximally 2.5% error
Cfg.opt_a.gref = 3;%This is reference G, for the initial optimization only. Ignore if not using newOpt

%NOTE: totalH is the time series for all subjects (nodes x timepoints)
%inside a 1 x numberOfSubjects array. sc is the C matrix.

subjID=num2str(getenv('subjID'));
load FC_CNPdata.mat
load SC_HCP_cortexOnly.mat
% Normalization that Gustavo always does to fit in with the global coupling value
C = AdjDens/max(AdjDens(:))*0.2;
subject_ts = squeeze(TS_ICA_AROMA_2P(:,:,subjID).');

[FC_simul, FC_emp, fitting, meta, ksP, PhasesD, Phases, bifpar]=...
    hbif_NewOpt(C, subject_ts, Tmax, wG, ldata, Cfg);

