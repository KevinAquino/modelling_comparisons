function saveeps(figH,fname,parameters,rend)

if(nargin<3)
    width = 4;height = 4;
else
    width = parameters(1);
    height = parameters(2);
end

if(nargin<4)
	figH.Renderer = 'painters';
else
	figH.Renderer = rend;
end

set(figH,'PaperPosition',[0.25 0.25 width height],'InvertHardCopy','off');
% set(figH,'Color',[0 0 0]);
% set(figH,'Color',[1 1 1]);
print(figH,[fname '.eps'],'-depsc');

