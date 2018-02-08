function res = poisson2d_Neumann(f)
globals;
    R                    = size(f,1); % nb lignes
    S                    = size(f,2); % nb colonnes

    R = Q+1; S = N+1;

    dn                   = 0:1:R-1;
    depn                 = 2*cos(pi*dn/R)-2;

    dm                   = 0:1:S-1;
    depm                 = 2*cos(pi*dm/S)-2;

    denom               = repmat(depn(:)*R^2,1,S) + repmat(depm(:)'*S^2,R,1);
    denom(denom(:)==0)  = 1;

    fhat                = mirt_dctn(f);
    uhat                = -(fhat)./denom;
    res                 = mirt_idctn(uhat);
end