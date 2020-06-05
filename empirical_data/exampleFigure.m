% 
noise_correction=3;
subject=84;



figure('color','white'); 
imagesc(corr(time_series(:,:,subject,noise_correction)))
set(gca,'YTick',[1:82],'YTickLabel',StructNames,'XTickLabelRotation',90)
set(gca,'XTick',[1:82],'XTickLabel',StructNames,'XTickLabelRotation',90)
axis image

title(['Subject: ',metadata.participant_id(subject,:),', Diagnosis: ',metadata.diagnosis(subject,:),'\newline age: ', ... 
	num2str(metadata.age(subject)),', meanFD: ',num2str(metadata.FD(subject)), ...
	', GSSD:',num2str(metadata.GSSD(subject)),', Noise Correction:', ...
	noiseOptions{noise_correction}],'fontSize',18)
 
 caxis([-0.5 0.5])
 colorbar



figure('color','white'); 
imagesc(((1:152)-1)*2,1:82,zscore(time_series(:,:,subject,noise_correction).',[],2))
set(gca,'YTick',[1:82],'YTickLabel',StructNames)
xlabel('time (s)')
colorbar;
caxis([-1.2,1.2]);colormap gray
title(['Subject: ',metadata.participant_id(subject,:),', Diagnosis: ',metadata.diagnosis(subject,:),'\newline age: ', ... 
	num2str(metadata.age(subject)),', meanFD: ',num2str(metadata.FD(subject)), ...
	', GSSD:',num2str(metadata.GSSD(subject)),', Noise Correction:', ...
	noiseOptions{noise_correction}],'fontSize',18)
% axis image


figure('color','white');
imagesc(ADJ_average)
set(gca,'YTick',[1:82],'YTickLabel',StructNames,'XTickLabelRotation',90)
set(gca,'XTick',[1:82],'XTickLabel',StructNames,'XTickLabelRotation',90)
axis image
title('ADJ Matrix','FontSize',18)
colorbar


% ROIS
figure('color','white')
plot3(roi_xyz(:,1),roi_xyz(:,2),roi_xyz(:,3),'.','MarkerSize',20)
title('ROI co-ordinates','FontSize',18)
axis off;