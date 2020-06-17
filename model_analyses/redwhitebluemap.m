% little function to make a colormap that goes from red to white to blue
function cmap=redwhitebluemap(Npoints)
	if(nargin<1)
		Npoints=100;
	end

	linearVector = linspace(0,1,Npoints).';
	reverseLinearVector = linearVector(end:-1:1);
	cmapA = [ones(size(linearVector)),reverseLinearVector,reverseLinearVector];
	cmapB = [linearVector,linearVector,ones(size(linearVector))];
	cmap = [cmapB;cmapA];