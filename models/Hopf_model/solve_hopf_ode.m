% This function is to make the Hopf simulation a little more universal, currently this line of code is repeated in some fashion.
function xs = solve_hopf_ode(omega,a,wC,dt,Tmax,TR,sig)
% omega: is the oscillatory frequency (in radians) for each node (calculated before hand)
% a    : is the hopf bifurcation parameter that transitions between oscillatory and noisy signals
% wC   : is the weighted coupling matrix
% dt   : is the time step (do we need some intuition behind this number? i guess it just has to converge)
% Tmax : is the max number of volumes in the fMRI time course
% TR   : is the Repetition Time for each volume
% sig  : is the strength of the noise terms -- currently at 0.1


% This here to work out the number of nodes
N = size(wC,1);

% Set up the degree matrix
sumC = repmat(sum(wC,2),1,2); % for sum Cij*xj

% Initialize the time series
xs=zeros(Tmax,N);

% Here solve the Hopf model using the parameter z which captures both x and y variables:
z = 0.1*ones(N,2); 
% z --> x = z(:,1), y = z(:,2)

% First perform the intial simulation, this is an SDE solved by the Euler–Maruyama method, where
% we have to include this step for solutions:
dsig = sqrt(dt)*sig;

% Need to check validity of the E-M method for this -- so far Deco et al. have used it for many papers but good to check with a different solver.


% Perform the first simulation, solving the SDE to initialize the simluation:
nn=0;
% discard first 3000 time steps
for t=0:dt:3000
    suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
    zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
    z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(N,2);
end

timeVector = 0:dt:((Tmax-1)*TR);
downsampleRate = TR/dt;
tsampled = (timeVector(1:downsampleRate:end));


% actual modeling (x=BOLD signal (Interpretation), y some other oscillation)
for tind=1:length(timeVector)
    suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
    zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
    z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(N,2);
    % Here we have the downsampling to the same TR -- i dont like this.
    if(ismember(timeVector(tind),tsampled))
        nn=nn+1;
        xs(nn,:)=z(:,1)';
    end
end