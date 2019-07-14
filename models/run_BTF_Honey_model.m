% This function runs the BTF model using the code straight from Honey and Sporns.

function [Vin,time,BOLD] = run_BTF_Honey_model(C,cparam)


% global variables (shared with 'simvec') (yuck)
global V1 V2 V3 V4 V5 V6 V6 V7 gCa gK gL VK VL VCa I b ani aei aie aee phi V8 V9 gNa VNa ane nse rnmda N CM vs c k_in

CM = C;
% runname
rn = ['modrun001_',num2str(cparam*100)];
init = 'randm';     % set if random initial condition is desired

N = size(CM,1);
% set out-strength (= out-degree for binary matrices)
k_in = sum(CM)';

% MODEL PARAMS =====================================
% set model parameters
V1 = -0.01; V2 = 0.15; V3 = 0; V4 = 0.3; V5 = 0; V7 = 0; V9 = 0.3; V8 = 0.15;
gCa = 1; gK = 2.0; gL = 0.5; gNa = 6.7;
VK = -0.7; VL = -0.5; I = 0.3; b = 0.1; phi = 0.7; VNa = 0.53; VCa = 1;
ani = 0.4; vs = 1; aei = 2; aie = 2; aee = 0.36; ane = 1; rnmda = 0.25;


%%%%%%ADJUST parameters
%phi = 0.2;
%%%%%%%%%%%%%%%%%%%%

% more parameters: noise, coupling, modulation
nse = 0;
c = cparam;            % ********* COUPLING ***********
modn = 0;
if (modn==0)
    V6 = 0.65;
else
    V6 = ones(N,1).*0.65 + modn*(rand(N,1)-0.5);
end;

% TIME PARAMS =====================================
% length of run and initial transient
% (in time segments, 1 tseg = l timesteps
tseg = 1;       % number of segments used in the intial transient
lseg = 10;      % number of segments used in the actual run
llen = 2000;   % length of each segment, in milliseconds
tres = 0.2;     % time resolution of model output, in milliseconds

% ERROR TOLERANCES =================================
% default: 'RelTol' 1e-3, 'AbsTol' 1e-6
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);

% INITIAL CONDITION ================================
% initial condition - random
if (strcmp(init,'randm'))
    ics = zeros(N,1);
    for i=1:N
        ics((i-1)*3+1:i*3,1) = [(rand-0.5)*0.8-0.2; (rand-0.5)*0.6+0.3; (rand-0.5)*0.16+0.05];
    end;
end;
% initial condition - start from an earlier run
if (strcmp(init,'saved'))
    load ics_previous    % substitute proper file name
end;

% START SIMULATION ================================
% TRANSIENT =======================================
disp('beginning dynamics (transient)...');
for SEGMENT=1:tseg
    tic;
    [t,y] = ode23('simvec',[0:tres:llen],ics,options);
    yics = y(end,:);
    disp(['finished segment ',num2str(SEGMENT)]);
    ics = yics;
    toc;
end;
disp('finished transient');

keyboard
% END TRANSIENT ==================================

% save model parameters and intial condition
eval(['save ',rn,'_params ics V1 V2 V3 V4 V5 V6 V6 V7 gCa gK gL VK VL VCa I b ani aei aie aee phi V8 V9 gNa VNa ane nse rnmda N CM vs c k_in']);

% RUN ============================================
% loop over 'lseg' segments of equal length
for SEGMENT=1:lseg
    tic;
    [t,y] = ode23('simvec',[0:tres:llen],ics,options);
    % keep only excitatory variable and downsample to 1 msec resolution
    Y = y(1:5:end,1:3:3*N-1);
    % save last time step as initial condition for next time segment
    yics = y(end,:);
    % save downsampled time series of excitatory variable, plus parameters
    eval(['save ',rn,'_part',num2str(SEGMENT),' Y yics V1 V2 V3 V4 V5 V6 V6 V7 gCa gK gL VK VL VCa I b ani aei aie aee phi V8 V9 gNa VNa ane nse rnmda N CM vs c k_in']);
    disp(['saved segment ',num2str(SEGMENT)]);
    % swap initial condition
    ics = yics;
    toc;
end;
% END OF RUN =======================================

% CONCATENATE OUTPUT FILES =========================
Yall = Yall_concat(rn,1:lseg);

% keyboard

Yall = Yall';
keyboard
% keyboard
for n=1:N,
    % get time series
    z = Yall(1:end,n);
    % Append the first 20 seconds to the front of the time series.
    ttrns = 20000;
    zt = Yall(1:ttrns+1,n);
    Vin(n,:) = [abs(diff([zt;z])); 0];

end



% Now comes the BOLD response..
dt = 1;
time = 0:dt:(size(Vin,2)-1)*dt;

% % Now model the HRF with a standard two-gamma distribution (can change this at any time..)
modelHrf = gampdf(time/1e3, 6, 1) - gampdf(time/1e3, 16, 1)/6;		
hrf = circshift(modelHrf.',round(length(modelHrf)/2)).';
% % Convolve the BOLD response with ms resolution
for n=1:N,
	BOLD(n,:) = conv(Vin(n,:),hrf,'same');
% 	BOLD_interp(n,:) = interp1(time(1e5:end)/1e3,zscore(BOLD(n,1e5:end)),TRs);
end

% After this we can interpolate etc..
