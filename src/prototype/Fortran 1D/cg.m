function x = cg(B,b,N,Q) 	
    x = zeros(size(b));
    r = b - B(x,N,Q);
    p = r;
    rold = r'*r;
    for i = 1:(Q+3)*(N+1)
        Ap = B(p,N,Q);
        alpha = rold/(p'*Ap);
        x = x + alpha*p;
        r = r - alpha*Ap;
        rnew = r'*r;
        if (sqrt(rnew) < 1e-10)  
            break 
        end 
        p = r + (rnew/rold)*p;
        rold = rnew;
    end 
end 