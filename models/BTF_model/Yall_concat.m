function  Yall = Yall_concat(rn,segs)

Yall = [];

% concatenate
for s=1:length(segs)
    eval(['load ',rn,'_part',num2str(segs(s))]);
    Yall = [Yall Y(1:end-1,:)'];
end;

% save Yall
eval(['save ',rn,'_Yall Yall']);

% delete the small segements...
for s=1:length(segs)
    eval(['delete ',rn,'_part',num2str(segs(s)),'.mat']);
end;
