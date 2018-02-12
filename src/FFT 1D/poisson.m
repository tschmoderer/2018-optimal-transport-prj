function p = poisson(f)
    globals;
    
    dn                    = 0:1:Q;
    depn                  = 2*cos(pi*dn/(Q+1))-2;

    dm                    = 0:1:N;
    depm                  = 2*cos(pi*dm/(N+1))-2;
    
    denom                = repmat(depn(:)*(Q+1)^2,1,N+1) + repmat(depm(:)'*(N+1)^2,Q+1,1);
    
    denom(denom(:)==0)  = 1;

    fhat                  = dct2(f);
    uhat                  = -(fhat)./denom;
    p                     = idct2(uhat);
end