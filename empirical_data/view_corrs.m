linearVector = linspace(0,1,100).';
reverseLinearVector = linearVector(end:-1:1);
cmapA = [ones(size(linearVector)),reverseLinearVector,reverseLinearVector];
cmapB = [linearVector,linearVector,ones(size(linearVector))];
cmap = [cmapB;cmapA];

for subject=1:120,
	close
	figure('color','white'); 
	subplot(1,3,1);
	imagesc(corr(time_series(:,:,subject,1)))
	title(['Subject: ',metadata.participant_id(subject,:), 'meanFD: ',num2str(metadata.FD(subject))],'fontSize',18)
	caxis([-1 1]);axis image;axis off
	subplot(1,3,2);
	imagesc(corr(time_series(:,:,subject,2)))
	caxis([-1 1]);axis image;axis off
	subplot(1,3,3);
	imagesc(corr(time_series(:,:,subject,3)))
	caxis([-1 1]);axis image;axis off
	colormap(cmap)
	pause
end
