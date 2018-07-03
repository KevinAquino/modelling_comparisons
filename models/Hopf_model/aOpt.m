 
function[aOptimized]=aOpt(a,Cfg,iN)

xs = zeros(3000/2,nNodes);
 wC = we*C;
 sumC = repmat(sum(wC,2),1,2); 
    
 
        minm=100;
        bestIter = 1; %JVS for tracking purposes
        for iter = 1:iN
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            z = 0.1*ones(nNodes,2); % --> x = z(:,1), y = z(:,2)
            nn=0;
            for t=1:dt:1000
                suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
                zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
                z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(nNodes,2);
            end
            
            for t=1:dt:3000%3000 seconds worth of simulated data; why dont we limit it to Tmax*Cfg.TRsec?
                suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
                zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
                z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(nNodes,2);
                
                if mod(t,Cfg.TRsec)==0
                    nn=nn+1;
                    xs(nn,:)=z(:,1)';
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            vsigs = zeros(1, nNodes);
            for seed=1:nNodes
                
                x=detrend(demean(xs(1:nn,seed)'));%was x
                ts_filt_wide =zscore(filtfilt(bfilt_wide,afilt_wide,x));%zscore before!
                
                TT=length(x);
                Ts = TT*Cfg.TRsec;
                freq = (0:TT/2-1)/Ts;
                [~, idxMinFreqS]=min(abs(freq-0.04));
                [~, idxMaxFreqS]=min(abs(freq-0.07));
                
                pw_filt_wide = abs(fft(ts_filt_wide));
                Pow1 = pw_filt_wide(1:floor(TT/2)).^2/(TT/2);
                Pow=gaussfilt(freq,Pow1,0.01);
                
                vsigs(seed)=sum(Pow(idxMinFreqS:idxMaxFreqS))/sum(Pow);
                
            end
            
            vsmin=min(vsigs);%vsig y vsigs max min points not equal
            vsmax=max(vsigs);
            bb=(vmax-vmin)/(vsmax-vsmin);
            aa=vmin-bb*vsmin;  %% adaptation of local bif parameters
            vsigs=aa+bb*vsigs;
            minm1=max(abs(vsig-vsigs)./vsig);
            trackminm1(iter, idx_g) = minm1;%JVS tracking of minm1
            
            if minm1<minm
                minm=minm1;
                a1=a;
                bestIter = iter; %JVS: tracking
                best_vsigs = vsigs; %JVS: tracking
            end
            
            %--------------------------------------------------------------
            %FEEDBACK
            %--------------------------------------------------------------
            fprintf(1, 'iter: %03d, minm1: %5.3f\n', iter, minm1);
            %--------------------------------------------------------------
            
            %CRITERION REACHED?
            if minm<Cfg.opt_a.abortCrit %default is 0.1
                break;
            end
            
            %UPDATE a VALUES FOR NEXT ITER
            if ~any(isnan(vsigs))
                a(:,1)=a(:,1)+Cfg.opt_a.updateStrength*(1-vsigs./vsig)';
                a(:,2)=a(:,2)+Cfg.opt_a.updateStrength*(1-vsigs./vsig)';
            else
                %JVS: this should not happen according to Gustavo, but it
                %can when a values run crazy
                warning('There are NaNs in the power spectra. Probably a-values too strong.');
                a = a1;
            end
        end
        aOptimized=a;
end