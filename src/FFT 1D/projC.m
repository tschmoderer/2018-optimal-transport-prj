function pU = projC(U)
    globals;
    
    U(1,1:N+1,2)   = f1;
    U(Q+2,1:N+1,2) = f0;

    p = (N+1)*diff(U(1:Q+1,:,1),[],2)  + (Q+1)*diff(U(:,1:N+1,2),[],1);

    f = poisson(-p);
    D = f;
    gf = -cat(3,(N+1)*[cat(2,-D(:,1),D(:,1:N) - D(:,2:N+1),D(:,N+1)); zeros(1,N+2)],(Q+1)*[cat(1,-D(1,:),D(1:Q,:) - D(2:Q+1,:),D(Q+1,:)),zeros(Q+2,1)]);
    
    pU = U; 
    pU(:,2:N+1,1) = pU(:,2:N+1,1) - gf(:,2:N+1,1);
    pU(2:Q+1,:,2) = pU(2:Q+1,:,2) - gf(2:Q+1,:,2);
end
