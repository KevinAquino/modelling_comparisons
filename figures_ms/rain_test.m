% figure;
% figure;
h1 = raincloud_plot(data2{1}, 'box_on', 1, 'color', theColors{1}, 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
     'box_col_match', 0);

h2 = raincloud_plot(data2{2}, 'box_on', 1, 'color', theColors{2}, 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);

h3 = raincloud_plot(data2{3}, 'box_on', 1, 'color', theColors{3}, 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .55, 'box_col_match', 0);
legend([h1{1} h2{1} h3{1}], noiseOptions);
% title(['Figure M7' newline 'A) Dodge Options Example 1']);
set(gca,'XLim', [-0.6 0.6], 'YLim', [-3 5],'fontSize',18,'YTick',[]);
% box off
axis tight

h=get(gca);
set(gca,'YLim',h.YLim + h.YLim/2);
% h.