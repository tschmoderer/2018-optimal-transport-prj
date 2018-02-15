function pU = projCs(U,V)
    globals; 
    
    b = U + interpAdj(V);
    pU = zeros(P+2,N+2,Q+2,3);
    
    r = b - pU - interpAdj(interp(pU));
    p = r;
    rold = sum(r(:).*r(:));
    
    for i = 1:2*(Q+2)*(N+2)*(P+2)
        Ip = p + interpAdj(interp(p));
        alpha = rold/sum(p(:).*Ip(:));
        pU = pU + alpha*p;
        r = r - alpha*Ip;
        rnew = sum(r(:).*r(:));
        if (sqrt(rnew) < 1e-10) 
            break 
        end
        p = r + (rnew/rold)*p;
        rold = rnew;
    end
end

