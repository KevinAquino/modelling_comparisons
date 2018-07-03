load hopf_scz.mat
bif_all_scz = bif_all;

load hopf_ctl

bif_all_parts = [bif_all bif_all_scz];


figure('Color','white');
imagesc(corr(bif_all_parts))
lbs = cell(10,1);lbs(1:5) = analyses;

for j=6:10,
	lbs{j} = ['SCZ_' analyses{j-5}];
end

set(gca,'XTickLabel',lbs,'YTickLabel',lbs,'TickLabelInterpreter','none','fontSize',18,'XTickLabelRotation',45);
caxis([0 1]);
title('Correlation of a')
colorbar
axis image
colormap jet