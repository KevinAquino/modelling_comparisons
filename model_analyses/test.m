IDScore=zeros(1,100);
for e=1:100;
    eucT1T2=zeros(200);euc=[];
    % Calculate the euclidean distance matrix. 
    for i=1:200;
        for j=1:200;
            euc=norm(d_wh_a(i,1:e)'-d_wh_b(j,1:e)'); 
            eucT1T2(i,j)=euc;
        end
    end
    % Calculate ID score per eigenvalue
    zeu=zscore(eucT1T2,0,'all');
    idx = find(~eye(size(zeu)));
    IDScore(e)=mean(zeu(idx))-mean(diag(zeu));
end