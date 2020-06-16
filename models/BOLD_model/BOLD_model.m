function bds = BOLD_model(z,time,TR)

	% time is in ms
	% TR is in s
	% z is the input into the convolution.


	% Here is the Boynton double gamma HRF
	modelHrf = gampdf(time/1e3, 6, 1) - gampdf(time/1e3, 16, 1)/6;		
	% Re-ordering the hrf so the convolution actually makes sense:
	hrf = circshift(modelHrf.',round(length(modelHrf)/2)).';

	% Here perform the convolutions.
	bds = conv(z,hrf,'same');

	% Now do an interpolation to get it all in the Time domain
	total_time=time(end)/1e3;
	new_time_vector=0:TR:total_time;

	% Now do simple linear interpolation (simple one)
	bds=interp1(time/1e3,bds,new_time_vector);

	