% Here define the general networks:
figure;
[square_mat,inds,total_order,all_regions] = nice_aparc_plotter(log10(C),[-4 -1],'white');

cols{4}=[0.470588,0.070588,0.513725];
cols{2}=[0.274510,0.509804,0.701961];
cols{5}=[0.000000,0.462745,0.050980];
cols{1}=[0.803922,0.243137,0.313725];
cols{3}=[0.862745,0.972549,0.647059];
cols{6}=[0.901961,0.580392,0.137255];
cols{7}=[0.768627,0.227451,0.984314];

cmap_aparc = zeros(82,3);
ii=1;
for nr=1:7,
	n_tot=length(inds{nr,1});
	factors=linspace(0.6,1,n_tot);
	for ind_r=1:n_tot
		cmap_aparc(ii,:) = factors(ind_r)*cols{nr};
		ii = ii+1;
	end

end

cmap_aparc(42:end,:) = cmap_aparc(1:41,:);

cmap_aparc = [0.75 0.75 0.75 ;cmap_aparc];

v2=1:82;
v2(total_order) = 1:82;


[left_figh,right_figh]=surface_display_results(label_struct,surf_left,surf_right,v2);
view([90 0])
colormap(cmap_aparc);caxis([0 82]);
set(left_figh,'visible','off');
right_figh.FaceColor='flat';
colorbar
[left_figh,right_figh]=surface_display_results(label_struct,surf_left,surf_right,v2);
view([-90 0])
colormap(cmap_aparc);caxis([0 82]);
set(left_figh,'visible','off')
right_figh.FaceColor='flat';
camlight
colorbar
