function pU = projCs(U,V)
    globals;
    b = U + interpAdj(V);
    pU = zeros(P+2,N+2,Q+1,3);
    r = b - pU - interpAdj(interp(pU));
    dir = r;
    rold = sum(r(:).*r(:));
    for i = 1:3*(P+2)*(N+2)*(Q+1)
        Idir = dir + interpAdj(interp(dir));
        alpha = rold/sum(dir(:).*Idir(:));
        pU = pU + alpha*dir;
        r = r - alpha*Idir;
        rnew = sum(r(:).*r(:));
        if (sqrt(rnew) < 1e-10) 
            break 
        end
        dir = r + (rnew/rold)*dir;
        rold = rnew;
    end
end

