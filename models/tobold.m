function out=tobold(zz,tt)
%  y = tobold(zz,tt)
% Solve BOLD ODEs for given input zz(tt) (tt in seconds), output bold y

% parameters
kappa = 0.65;   % rate of signal decay, s^-1                0.65
gamma = 0.41;   % rate of flow-dependent elimination, s^-1  0.41
tau = 0.98;     % haemodynamic transit time, s              0.98
alpha = 0.32;   % Grubb's exponent                          0.32
rho = 0.34;     % resting oxygen extraction fraction        0.34
V0 = 0.02;      % resting blood volume fraction             0.02


% ICs
%ics = [-0.0004    3.4415    1.4851    0.4968]'; % ?
%ics=[0 1 1 1];
%z0=zz(1);
z0=mean(zz);
f0=z0/gamma+1;
v0=f0^alpha;
q0=f0*(1-(1-rho)^(1/f0))/(rho*v0^(1/alpha-1));
ics=[0;f0;v0;q0];

% interpolant
F = griddedInterpolant(tt,zz);

% solve
tzero = tt(1);
tend = tt(end);
%tsteps = 10*tend;
%[t,x] = rgk4('boldodes_yy',tzero,tend,ics,(tend-tzero)*tsteps);
sol = ode45(@boldodes_yy,[tzero tend],ics,odeset('MaxStep',0.002));
bolddt=0.1;
t=tzero:bolddt:tend;
x=deval(sol,t)';

% BOLD signal
k1 = 7*rho;
k2 = 2;
k3 = 2*rho - 0.2;
y = V0*(k1*(1-x(:,4))+k2*(1-x(:,4)./x(:,3))+k3*(1-x(:,3)));

out.y=y;
out.t=t;
% keyboard

function fn=boldodes_yy(t,x)
% Haemodynamic model embedding the Balloon-Windkessel model
% allowing input time series z(tt)
% x(1) = s = vasodilatory signal
% x(2) = f = inflow
% x(3) = v = blood volume
% x(4) = q = deoxyhaemoglobin content
% z = neuronal activity

%z=interp1(tt,zz,t);
z=F(t);
fn=[z-kappa*x(1) - gamma*(x(2)-1);
    x(1);
    (x(2) - (x(3)).^(1/alpha))/tau;
    (x(2).*(1-(1-rho).^(1/x(2)))/rho - (x(3)).^((1-alpha)/alpha).*x(4))/tau;];
end




end
