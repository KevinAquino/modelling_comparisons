% testing

t_max = 1e4;
dt = 1e-2;
t = 0:dt:t_max;

beta = 0.02;
prevx = 1;
prevy = 1;
f = 0.006;
omega = 2*pi*f;
val = 0;
% Hopf parameter
for a = linspace(-1,1,20);
	disp(a);
	val = 1+val;
	% Newtonian method here for the Hopf Model solution -- test
	for nt=2:length(t)
		x(nt) = prevx + dt*prevx*(a-prevx^2 - prevy^2) - omega*prevy*dt + sqrt(dt)*beta*randn(1);
		y(nt) = prevy + dt*prevy*(a-prevx^2 - prevy^2) + omega*prevx*dt + sqrt(dt)*beta*randn(1);
		prevx = x(nt);
		prevy = y(nt);
	end	
	x_all(:,val) = x;
end
