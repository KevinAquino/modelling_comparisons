%Function to create linear interpolation of bifurcation parameters. Input
%arguments are the reference a parameter vector, reference coupling G and
%the coupling vector G. Output is a Matrix bpN, which rows contain
%linearized parameter values and columns the number of nodes.
%Victor Saenger, 2016.


function[bpN]=aLin(a,gr,wG)

%Normalize values between -1 and 1:
gpN=zeros(1,length(a));
pi=find(a>0);
ni=find(a<0);

pgp=a(pi);
ngp=a(ni);
n_pgp=pgp/max(pgp);
n_ngp=ngp/abs(min(ngp));

gpN(pi)=n_pgp;
gpN(ni)=n_ngp; 

% Scale parameter assuming linear evolution:
p=(min(wG)/gr):(abs(wG(1)-wG(2))/gr):(max(wG)/gr);%this is the new sclaing vector
bpN=zeros(length(p),length(a));

    for n = 1:length(gpN);
        if gpN(n)<0
            bpN(:,n)=gpN(n)*p*abs(min(a));
        else
            bpN(:,n)=gpN(n)*p*max(a);
        end
    end

end
