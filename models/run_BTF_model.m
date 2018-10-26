% This here runs the BTF model while making sure we have the right parameters. All the parameters are default within the BDT
% albeit for a handful of very important parameters! These made all the difference actually, so here we are making them
% explicitly defined in this regime. 

function [Vin,BOLD] = run_BTF_model(C,cparam)
BOLD=[];

% Here is how long each simulation segment is being run for, there is no point doing the whole thing as it will kill memory
simSegmentTime = 2000;

sys = BTF2003(C);
sys.pardef=bdSetValue(sys.pardef,'C',cparam);

deltaV = bdGetValue(sys.pardef,'deltaV');
deltaZ = bdGetValue(sys.pardef,'deltaZ');


% Change these parameters because these are the ones used by Honey et al. 
sys.pardef=bdSetValue(sys.pardef,'deltaV',0.65*ones(size(deltaV)));
sys.pardef=bdSetValue(sys.pardef,'deltaZ',0.65*ones(size(deltaZ)));

% Change the initial conditions to be consistent with prev. work
nnodes=size(C,1);
randvec=rand(nnodes,1);
sys.vardef = [ struct('name','V', 'value',(randvec-0.5)*0.8-0.2);      % Mean firing rate of excitatory cells
               struct('name','W', 'value',(randvec-0.5)*0.6+0.3);      % Proportion of open K channels
               struct('name','Z', 'value',(randvec-0.5)*0.16+0.05) ];  % Mean firing rate of inhibitory cells

Gion = bdGetValue(sys.pardef,'Gion');
Gion(1) = 1; %GCa has to be this value for Honey et al. paper
sys.pardef = bdSetValue(sys.pardef,'Gion',Gion);

aee = 0.36;% This here is the value from honey et al. 
sys.pardef = bdSetValue(sys.pardef,'aee',aee);

% Solve the model for 2 seconds to remove transients
sol = bdSolve(sys,[0 simSegmentTime],@ode23);

Y = reshape(sol.y(:,end),size(C,1),3,1);
for j=1:3,	
	sys.vardef(j).value = Y(:,j);
end

% After this, let the system run.
sol = bdSolve(sys,[0 simSegmentTime],@ode23);
for j=1:size(C,1)*3,
	y_all(j,:) = interp1(sol.x,sol.y(j,:),[1:1:simSegmentTime]);
end

Y = reshape(y_all,size(C,1),3,simSegmentTime);

% Take the last values
lastValues = Y(:,:,end);
Vtest = squeeze(Y(:,1,:));

% Take the deriviate of the excitatory potentials (this is what has been used in the past)
for n=1:size(C,1),	
	inputNeural(n,:) = abs(gradient(Vtest(n,:),1));	
end

Vin = inputNeural;


Vtin = [];


% Solve this for 200 more seconds;
% There is a tool within BDToolbox that does this, but here I want more control to 
% observe all the parameters and not save all the values.

for n=1:200,	
	disp(['Iteration: ',num2str(n),'....']);
	% Now set the last point as the intial value for the simulation to continue the simulation
	sys2 = sys;
	for j=1:3,
		sys2.vardef(j).value = lastValues(:,j);
	end
	% Solve again for 2 more seconds
	sol2 = bdSolve(sys2,[0 simSegmentTime],@ode23);

	for j=1:size(C,1)*3,
		y_all(j,:) = interp1(sol2.x,sol2.y(j,:),[1:1:simSegmentTime]);
	end


	Y = reshape(y_all,size(C,1),3,simSegmentTime);


	Vtest = squeeze(Y(:,1,:));

	for n=1:size(C,1),		
		inputNeural(n,:) = abs(gradient(Vtest(n,:),1));		
	end

	Vin = [Vin,inputNeural];
	Vtin = [Vtin,Vtest];
	lastValues = Y(:,:,end);

end


% Downsample the responses to 0.1s resolution for the responses and peform a convolution. Note that the Balloon-model
% is not really needed, as the responses are in the linear regime, it just adds complexity to the modelling, and although
% the Balloon model displays some comparable characteristics of the BOLD response, it is not a model of reality. As shown
% by Hillman et al 2007, and should be replaced by a better model e.g. Aquino et al. 2012, Ress et al 2013.
dt = 1;

% I have downsampled the responses here rather crudely -- need to do this properly by filtering then downsampling.
% not sure if this will make a difference. 
dsRate=100;

time = 0:dsRate*dt:(size(Vin,2)-1)*dt;

% Here is the Boynton double gamma HRF
modelHrf = gampdf(time/1e3, 6, 1) - gampdf(time/1e3, 16, 1)/6;		
% Re-ordering the hrf so the convolution actually makes sense:
hrf = circshift(modelHrf.',round(length(modelHrf)/2)).';

% Here perform the convolutions.
for n=1:size(C,1),
	BOLD(n,:) = conv(Vin(n,1:dsRate:end),hrf,'same');
end


