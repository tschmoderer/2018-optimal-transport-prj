function [pmC,pfC] = projC(mbar,fbar,N,Q)
    normalise = @(f) f/sum(f(:));

    f0 = normalise(gauss(0.1,0.005,N));
    f1 = normalise(gauss(0.9,0.005,N));
    
    Cst = poisson(zeros(Q+1,N+1),0,0,f1,f0,N,Q,1e-5);
    Cstmbar = zeros(Q+1,N+2); Cstfbar = zeros(Q+2,N+1);
    
    Cstmbar(:,1) = -Cst(:,1); Cstmbar(:,end) = Cst(:,end);
    Cstmbar(:,2:end-1) = Cst(:,1:end-1) - Cst(:,2:end);

    Cstfbar(1,:) = -Cst(1,:); Cstfbar(end,:) = Cst(end,:); 
    Cstfbar(2:end-1,:) = Cst(1:end-1,:) - Cst(2:end,:);
    
    Cstmbar = N*Cstmbar; Cstfbar = Q*Cstfbar;
    
    d = N*(mbar(:,2:end) - mbar(:,1:end-1)) + Q*(fbar(2:end,:) - fbar(1:end-1,:));
    
   % S = poisson(d,d(:,1),d(:,end),d(1,:),d(end,:),N,Q,1e-5);
    S = poisson(d,mbar(:,1),mbar(:,end),fbar(1,:),fbar(end,:),N,Q,1e-5);
    
    Smbar = zeros(Q+1,N+2); Sfbar = zeros(Q+2,N+1);
    
    Smbar(:,1) = -S(:,1); Smbar(:,end) = S(:,end);
    Smbar(:,2:end-1) = S(:,1:end-1) - S(:,2:end);
    
    Sfbar(1,:) = -S(1,:); Sfbar(end,:) = S(end,:); 
    Sfbar(2:end-1,:) = S(1:end-1,:) - S(2:end,:);
    
    Smbar = N*Smbar; Sfbar = Q*Sfbar;
    
    pmC = mbar - Smbar + Cstmbar;
    pfC = fbar - Sfbar + Cstfbar;
end

