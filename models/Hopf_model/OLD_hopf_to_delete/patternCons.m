
function [phFCD meta meanp stdp]=patternCons(Phases,N,Tmax)

Isubdiag = find(tril(ones(N),-1));

 %phFCD=zeros(N,Tmax);

     T=10:Tmax-10;
     
     
     %pattern=zeros(length(T-9),N);
     for t=T
      kudata=sum(complex(cos(Phases(:,t)),sin(Phases(:,t))))/N;
      syncdata(t-9)=abs(kudata);
      for i=1:N
        for j=1:i-1
         patt(i,j)=cos(adif(Phases(i,t),Phases(j,t)));
        end
      end 
      pattern(t-9,:)=patt(Isubdiag);
     end
     
      meta=std(syncdata);

     kk3=1;
     npattmax=size(pattern,1);
     for t=1:npattmax-2
       p1=mean(pattern(t:t+2,:));
       %phfcddata=zeros(1,N);
      for t2=t+1:npattmax-2
       p2=mean(pattern(t2:t2+2,:));
       phFCD(kk3)=dot(p1,p2)/norm(p1)/norm(p2);
       kk3=kk3+1;
      end
     end 
     %metastabilitydata2(nsub)=std(syncdata);


    meanp=mean(phFCD);
    stdp=std(phFCD);
end

