% Left Hemi
[vertices, label, ctab] = read_annotation('lh.HCP-MMP1.annot');
[verticesr, labelr, ctabr] = read_annotation('rh.HCP-MMP1.annot');

[vl,~] = read_surf('/Applications/freesurfer/subjects/fsaverage/surf/lh.white');
[vr,~] = read_surf('/Applications/freesurfer/subjects/fsaverage/surf/rh.white');

for j=1:180;
	parc1 = vl(find(label==ctab.table(j+1,5)),:);
	vertexParcel(j,:) = mean(parc1,1);
end

for j=1:180;
	parc1 = vr(find(labelr==ctabr.table(j+1,5)),:);
	vertexParcel(j+180,:) = mean(parc1,1);
end

% Distance Adj Matrix
for j=1:360,
	for nc=j:360
		distMatrix(j,nc) = sqrt(sum(abs(vertexParcel(j,:) - vertexParcel(nc,:)).^2));
	end
end

distMatrix = distMatrix + distMatrix.' - eye(size(distMatrix));

% Here make the connection stuff
factor = 0.4;
Adj = zeros(size(distMatrix));
for j=1:360,
	for nc=j:360
		hh = rand;
		prob = exp(-factor*distMatrix(j,nc)/0.01e3);
		if(hh<prob)
			Adj(j,nc)=prob;
		end
	end
end

Adj = Adj + Adj.' - eye(size(Adj));
figure;imagesc(Adj)