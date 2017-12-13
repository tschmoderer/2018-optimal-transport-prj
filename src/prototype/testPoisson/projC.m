function [pmC,pfC] = projC(mbar,fbar)
    globals;

    %Cst = poisson(zeros(Q+1,N+1),0,0,f1,f0,N,Q,1e-5);
    % Calcul de la constante %
    %!!%
    ymbar = zeros(Q+1,N+2); 
    yfbar = zeros(Q+2,N+1); yfbar(1,:) = f1; yfbar(end,:) = f0;
    yy = N*(ymbar(:,2:end)-ymbar(:,1:end-1)) + Q*(yfbar(2:end,:) - yfbar(1:end-1,:));
    %!!%

    C1 = poisson2d_Neumann(-yy); 
    Cmbar = zeros(Q+1,N+2); Cfbar = zeros(Q+2,N+1);

    Cmbar(:,1) = -C1(:,1); Cmbar(:,end) = C1(:,end);  Cmbar(:,1) = 0; Cmbar(:,end) = 0;
    Cmbar(:,2:end-1) = (C1(:,1:end-1) - C1(:,2:end));

    Cfbar(1,:) = -C1(1,:); Cfbar(end,:) = C1(end,:); Cfbar(1,:) = f1/Q; Cfbar(end,:) = f0/Q;
    Cfbar(2:end-1,:) = (C1(1:end-1,:) - C1(2:end,:));

    Cstmbar = N*Cmbar; Cstfbar = Q*Cfbar; 
    
    % calcul de la projection
    
    d = N*(mbar(:,2:end) - mbar(:,1:end-1)) + Q*(fbar(2:end,:) - fbar(1:end-1,:));

    S = poisson2d_Neumann(-d);

    Smbar = zeros(Q+1,N+2); Sfbar = zeros(Q+2,N+1);
    
    Smbar(:,1) = mbar(:,1)/N; Smbar(:,end) = mbar(:,end)/N;
    Smbar(:,2:end-1) = S(:,1:end-1) - S(:,2:end);
    
    Sfbar(1,:) = fbar(1,:)/Q; Sfbar(end,:) = fbar(end,:)/Q;
    Sfbar(2:end-1,:) = S(1:end-1,:) - S(2:end,:);
    
    Smbar = N*Smbar; Sfbar = Q*Sfbar;

    % la projection
    pmC = mbar - Smbar + Cstmbar;
    pfC = fbar - Sfbar + Cstfbar;
end

