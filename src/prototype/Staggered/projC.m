function pU = projC(U)
    globals;
    y = zeros(P+3,N+3,Q+3);
    y(1:P+1,1:N+1,Q+2) = f0;
    y(1:P+1,1:N+1,Q+3) = f1;

    x = zeros(P+3,N+3,Q+3);
    b = y - A(U);
    r = b - A(AS(x));
    p = r;
    rold = sum(r(:).*r(:));
    for i = 1:(Q+3)*(N+3)*(P+3)
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
       

%       
%   U(1:P+1,1:N+1,1,3)   = f0;
%   U(1:P+1,1:N+1,Q+2,3) = f1;
% 
%   U([1 end],:,:,1) = 0; U(:,[1 end],:,2) = 0;
% 
%    D = (P+1)*diff(U(:,1:N+1,1:Q+1,1),[],1) + (N+1)*diff(U(1:P+1,:,1:Q+1,2),[],2) + (Q+1)*diff(U(1:P+1,1:N+1,:,3),[],3);
%    p = poisson(-D);
%    pU = U;
% 
%   pU(2:(end-1),1:N+1,1:Q+1,1)   = U(2:(end-1),1:N+1,1:Q+1,1) - diff(p,[],1)*(P+1);
%   pU(1:P+1,2:(end-1),1:Q+1,2)   = U(1:P+1,2:(end-1),1:Q+1,2) - diff(p,[],2)*(N+1);
%   pU(1:P+1,1:N+1,2:(end-1),3)   = U(1:P+1,1:N+1,2:(end-1),3) - diff(p,[],3)*(Q+1);
end

