function pU = projC(U)
    globals;
   
% 
%     x = zeros(N+3,P+3,Q+3);
%     b = y - A(U);
%     r = b - A(AS(x));
%     p = r;
%     rold = sum(r(:).*r(:));
%     for i = 1:(N+3)*(P+3)*(Q+3)
%         Ap = A(AS(p));
%         alpha = rold/sum(p(:).*Ap(:));
%         x = x + alpha*p;
%         r = r - alpha*Ap;
%         rnew = sum(r(:).*r(:));
%         if (sqrt(rnew) < 1e-10) 
%             break
%         end
%         p = r + (rnew/rold) * p;
%         rold = rnew;
%     end
% 
%     pU = U + AS(x);
       


    U(1:N+1,1:P+1,1  ,3) = f1;
    U(1:N+1,1:P+1,Q+2,3) = f0;
    U([1 end],:,:,1) = 0; U(:,[1 end],:,2) = 0;

    D = (N+1)*diff(U(:,1:P+1,1:Q+1,1),[],1) + (P+1)*diff(U(1:N+1,:,1:Q+1,2),[],2) + (Q+1)*diff(U(1:N+1,1:P+1,:,3),[],3);

    p = poisson(-D);


    pU2 = U; 

    pU2(2:N+1,1:P+1,1:Q+1,1) = U(2:N+1,1:P+1,1:Q+1,1) - (N+1)*diff(p,[],1);
    pU2(1:N+1,2:P+1,1:Q+1,2) = U(1:N+1,2:P+1,1:Q+1,2) - (P+1)*diff(p,[],2);
    pU2(1:N+1,1:P+1,2:Q+1,3) = U(1:N+1,1:P+1,2:Q+1,3) - (Q+1)*diff(p,[],3);


%     [sum(sum(sum(abs(D-div(U))))) sum(sum(sum(sum(abs(pU-pU2)))))]


pU = pU2;
end

