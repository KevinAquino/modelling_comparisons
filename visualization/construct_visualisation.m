[vertices,faces]=read_surf('lh.inflated');
vertices(:,1) = -45+vertices(:,1);
surf_left.vertices=vertices;
surf_left.faces=faces+1;
surf_left.Edgecolor='none';

[vertices,faces]=read_surf('rh.inflated');
vertices(:,1) = 45+vertices(:,1);
surf_right.Edgecolor='none';
surf_right.faces=faces+1;
surf_right.vertices=vertices;


[~, label, ctab] = read_annotation('lh.aparc.annot');

label_struct.left_label=label;
label_struct.left_ctab=ctab;

[~, label, ctab] = read_annotation('rh.aparc.annot');
label_struct.right_label=label;
label_struct.right_ctab=ctab;