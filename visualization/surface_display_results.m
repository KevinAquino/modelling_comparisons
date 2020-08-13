
% A function to display results on a surface
function [left_figh,right_figh]=surface_display_results(label_struct,surf_left,surf_right,vector)
if(nargin<4)
	vector=[1:82];
end

if(nargin<5)
	view_angle=[45 45];
end
% MAKE NANs zero
vector(isnan(vector))=0;

% Paint each surface with the vector colour
% LEFT
parcel_indices=1+setdiff([1:35],4);
surf_left.CData=zeros(size(surf_left.vertices,1),1);
for ind=1:length(parcel_indices),
	parcel_vertex_inds=find(label_struct.left_label==label_struct.left_ctab.table(parcel_indices(ind),5));
	surf_left.CData(parcel_vertex_inds) = vector(ind);
end

surf_right.CData=zeros(size(surf_right.vertices,1),1);
for ind=1:length(parcel_indices),
	parcel_vertex_inds=find(label_struct.right_label==label_struct.right_ctab.table(parcel_indices(ind),5));
	surf_right.CData(parcel_vertex_inds) = vector(41+ind);
end

figure('color','white');hold on;
left_figh=patch(surf_left,'FaceColor','interp');
right_figh=patch(surf_right,'FaceColor','interp');
view(view_angle);

% Special patch properties
camlight
material dull
axis image
axis off;