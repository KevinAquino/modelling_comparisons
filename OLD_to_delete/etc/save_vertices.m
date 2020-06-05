% Save vertices_matlab

[vertices,~] = read_surf(['/Applications/freesurfer/subjects/fsaverage/surf/lh.white']);
[vertex,label,ctab] = read_annotation('/Applications/freesurfer/subjects/fsaverage/label/lh.aparc.annot');

for j=2:length(ctab.struct_names)
	inds = find(label==ctab.table(j,5));
	position_lh(j-1,:) = mean(vertices(inds,:),1);
end

[vertices,~] = read_surf(['/Applications/freesurfer/subjects/fsaverage/surf/rh.white']);
[vertex,label,ctab] = read_annotation('/Applications/freesurfer/subjects/fsaverage/label/rh.aparc.annot');

for j=2:length(ctab.struct_names)
	inds = find(label==ctab.table(j,5));
	position_rh(j-1,:) = mean(vertices(inds,:),1);
end

centroids = [position_lh;position_rh];

centroids = centroids([1:3,5:38,40:end],:);

figure;plot3(centroids(:,1),centroids(:,2),centroids(:,3),'.','MarkerSize',20)




% Here now it is for HCP

[vertices,~] = read_surf(['/Applications/freesurfer/subjects/fsaverage/surf/lh.white']);
[vertex,label,ctab] = read_annotation('/Applications/freesurfer/subjects/fsaverage/label/lh.hcp-MMP1.annot');

for j=2:length(ctab.struct_names)
	inds = find(label==ctab.table(j,5));
	position_lh(j-1,:) = mean(vertices(inds,:),1);
end

[vertices,~] = read_surf(['/Applications/freesurfer/subjects/fsaverage/surf/rh.white']);
[vertex,label,ctab] = read_annotation('/Applications/freesurfer/subjects/fsaverage/label/rh.hcp-MMP1.annot');

for j=2:length(ctab.struct_names)
	inds = find(label==ctab.table(j,5));
	position_rh(j-1,:) = mean(vertices(inds,:),1);
end


centroids = [position_lh;position_rh];

load('scz_ten.mat');

save('scz_hcp_ten','time_series_hcp','centroids','analyses');