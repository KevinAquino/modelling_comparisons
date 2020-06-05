function Joptim=Balance_J(we,C)

  % Parameters for the simulations
  Nnew=size(C,1);
  dt=0.1;
  tmax1=10000;
  tspan1=0:dt:tmax1;
  % Load up all the parameters in a seperate file.  
  load_BEI_parameters;
  % Initial estimates for J
  J=Balance_J_analytic(we,C);
  % Test  
  curr=zeros(tmax1,Nnew);
  % Change in the parameters:
  delta=delta_J*ones(Nnew,1);

  for k=1:50000
   % Two initialized conditions
   sn=0.001*ones(Nnew,1);
   sg=0.001*ones(Nnew,1);
   nn=1;
   j=0;
   for i=2:1:length(tspan1)
    xn=I0*Jexte+w*JN*sn+we*JN*C*sn-J.*sg;
    xg=I0*Jexti+JN*sn-sg;
    rn=phie(xn);
    rg=phii(xg);
    sn=sn+dt*(-sn/taon+(1-sn)*gamma_factor.*rn./1000.)+sqrt(dt)*sigma_factor*randn(Nnew,1);
    sn(sn>1) = 1;  
    sn(sn<0) = 0;      
    sg=sg+dt*(-sg/taog+rg./1000.)+sqrt(dt)*sigma_factor*randn(Nnew,1);
    sg(sg>1) = 1;        
    sg(sg<0) = 0;
    j=j+1;
    if j==10
     curr(nn,:)=xn'-I_phie/c_phie;
     nn=nn+1;
     j=0;
    end
   end

   currm=mean(curr(1000:end,:),1);
  % keyboard
   flag=0;
   for n=1:1:Nnew
    if abs(currm(n)+I_condition)>tolerance_I
     if currm(n)<-I_condition 
      J(n)=J(n)-delta(n);
      % Optimzation limits
      delta(n)=delta(n)-optimization_limit_J;
      if delta(n)<optimization_limit_J;
         delta(n)=optimization_limit_J;
      end
     else 
      J(n)=J(n)+delta(n);
     end
    else
     flag=flag+1;
    end
   end

   if flag==Nnew
    break;
   end
  end

  Joptim=J;
  % keyboard
  end
