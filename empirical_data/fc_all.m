 
for k=1:3,
	for j=1:120,
	 	FC(:,:,j,k) = corr(time_series(:,:,j,k));%                152x82x260x3
	 end
end



figure('color','white');
subject=15;

subplot(131)
imagesc(FC(:,:,subject,1));axis image;colormap(cmap);axis off;caxis([-1 1])
subplot(132)
imagesc(FC(:,:,subject,2));axis image;colormap(cmap);axis off;caxis([-1 1])
subplot(133)
imagesc(FC(:,:,subject,3));axis image;colormap(cmap);axis off;caxis([-1 1])