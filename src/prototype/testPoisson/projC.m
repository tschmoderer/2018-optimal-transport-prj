function [pmC,pfC] = projC(mbar,fbar,N,Q)
    normalise = @(f) f/sum(f(:));

    f0 = normalise(0.01 + gauss(0.4,0.05,N));
    f1 = normalise(0.01 + gauss(0.6,0.05,N));
    
    %Cst = poisson(zeros(Q+1,N+1),0,0,f1,f0,N,Q,1e-5);
    y = zeros(Q+1,N+1); y(1,:) = Q*f1; y(end,:) = -Q*f0; %!!%
    Cst = poisson2d_Neumann(y);
    Cstmbar = zeros(Q+1,N+2); Cstfbar = zeros(Q+2,N+1);
    
    Cstmbar(:,1) = -Cst(:,1); Cstmbar(:,end) = Cst(:,end); Cstmbar(:,1) = 0; Cstmbar(:,end) = 0; %!!% 
    Cstmbar(:,2:end-1) = Cst(:,1:end-1) - Cst(:,2:end);

    Cstfbar(1,:) = -Cst(1,:); Cstfbar(end,:) = Cst(end,:);  Cstfbar(1,:) = f1/Q; Cstfbar(end,:) = f0/Q; %!!%
    Cstfbar(2:end-1,:) = Cst(1:end-1,:) - Cst(2:end,:);
    
    Cstmbar = N*Cstmbar; Cstfbar = Q*Cstfbar;
    
    d = N*(mbar(:,2:end) - mbar(:,1:end-1)) + Q*(fbar(2:end,:) - fbar(1:end-1,:));
    
   % S = poisson(d,d(:,1),d(:,end),d(1,:),d(end,:),N,Q,1e-5);
   % S = poisson(d,mbar(:,1),mbar(:,end),fbar(1,:),fbar(end,:),N,Q,1e-5);
    %!!%
    d(1,:) = Q*fbar(1,:); d(end,:) = -Q*fbar(end,:);
    d(:,1) = N*mbar(:,1); d(:,end) = -N*mbar(:,end);
    S = poisson2d_Neumann(-d);
    %!!%
    Smbar = zeros(Q+1,N+2); Sfbar = zeros(Q+2,N+1);
    
    Smbar(:,1) = -S(:,1); Smbar(:,end) = S(:,end); Smbar(:,1) = mbar(:,1)/N; Smbar(:,end) = mbar(:,end)/N; %!!%
    Smbar(:,2:end-1) = S(:,1:end-1) - S(:,2:end);
    
    Sfbar(1,:) = -S(1,:); Sfbar(end,:) = S(end,:); Sfbar(1,:) = fbar(1,:)/Q; Sfbar(end,:) = fbar(end,:)/Q; %!!%
    Sfbar(2:end-1,:) = S(1:end-1,:) - S(2:end,:);
    
    Smbar = N*Smbar; Sfbar = Q*Sfbar;
%     %!!% 
%     Smbar(:,2:end-1) = Smbar(:,2:end-1)/N;
%     Sfbar(2:end-1,:) = Sfbar(2:end-1,:)/Q;
%     %!!!%
%     
%     Smbar = N*Smbar; Sfbar = Q*Sfbar;

    pmC = mbar - Smbar + Cstmbar;
    pfC = fbar - Sfbar + Cstfbar;
end

