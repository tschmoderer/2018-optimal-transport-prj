function p = poisson(f)
%    globals;
    
    N                     = size(f,1);
    M                     = size(f,2);
    R                     = size(f,3);
    hx                    = 1/N;
    hy                    = 1/M;
    hz                    = 1/R;

    dn                    = 0:1:N-1;
    depn                  = 2*cos(pi*dn/N)-2;
    dm                    = 0:1:M-1;
    depm                  = 2*cos(pi*dm/M)-2;
    dr                    = 0:1:R-1;
    depr                  = 2*cos(pi*dr/R)-2;

    denom                = repmat(depn(:)/hx^2,[1,M,R]) + repmat(depm(:)'/hy^2,[N,1,R]) + permute(repmat(depr(:)/hz^2,[1,N,M]),[2 3 1]);
    denom(denom(:)==0)  = 1;



    fhat                  = dctn(f);
    uhat                  = -(fhat)./denom;
    p                     = idctn(uhat);
end