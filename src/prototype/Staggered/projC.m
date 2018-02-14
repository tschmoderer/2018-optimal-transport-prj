function pU = projC(U)
    globals;
    y = zeros(Q+3,N+3);
    y(Q+2,1:N+1) = f0;
    y(Q+3,1:N+1) = f1;

    x = zeros(Q+3,N+3);
    b = y - A(U);
    r = b - A(AS(x));
    p = r;
    rold = sum(r(:).*r(:));
    for i = 1:(Q+3)*(N+3)
        Ap = A(AS(p));
        alpha = rold/sum(p(:).*Ap(:));
        x = x + alpha*p;
        r = r - alpha*Ap;
        rnew = sum(r(:).*r(:));
        if (sqrt(rnew) < 1e-10) 
            break
        end
        p = r + (rnew/rold) * p;
        rold = rnew;
    end

    pU = U + AS(x);
end

