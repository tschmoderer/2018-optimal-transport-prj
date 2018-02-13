function pU = projC(U)
    globals;
    y = zeros(P+3,N+3,Q+2);
    y(1:P+1,1:N+1,Q+1) = f0;
    y(1:P+1,1:N+1,Q+2) = f1;
    x = zeros(P+3,N+3,Q+2);
    b = y - A(U);
    r = b - A(AS(x));
    dir = r;
    rold = sum(r(:).*r(:));
    for i = 1:(Q+3)*(N+3)
        Adir = A(AS(dir));
        alpha = rold/sum(dir(:).*Adir(:));
        x = x + alpha*dir;
        r = r - alpha*Adir;
        rnew = sum(r(:).*r(:));
        if (sqrt(rnew) < 1e-10) 
            break
        end
        dir = r + (rnew/rold) * dir;
        rold = rnew;
    end
		
    pU = U + AS(x);
    Apu = A(pU);
    err = norm(y(:) - Apu(:))/norm(y(:))
    %% gniah 
    
%     U(1:P+1,1:N+1,1,3)   = f1;  U(1:P+1,1,1:Q,1)   = 0; U(1,1:N+1,1:Q,2)   = 0;
%     U(1:P+1,1:N+1,Q+1,3) = f0;  U(1:P+1,N+2,1:Q,1) = 0; U(P+2,1:N+1,1:Q,2) = 0;
% 
%     p = (N+1)*diff(U(1:P+1,:,1:Q,1),[],2)  + (P+1)*diff(U(:,1:N+1,1:Q,2),[],1)  + (Q)*diff(U(1:P+1,1:N+1,:,3),[],3);
% 
%     f = poisson(-p);
%     
%     gf = zeros(P+2,N+2,Q+1,3);
%     gf(1:P+1,1,1:Q,1)     = -f(:,1,:);
%     gf(1:P+1,2:N+1,1:Q,1) = f(:,1:N,:) - f(:,2:N+1,:);
%     gf(1:P+1,N+2,1:Q,1)   = f(:,N+1,:);
%     gf(:,:,:,1) = (N+1)*gf(:,:,:,1);
% 
%     gf(1,1:N+1,1:Q,2)     = -f(1,:,:);
%     gf(2:P+1,1:N+1,1:Q,2) = f(1:P,:,:) - f(2:P+1,:,:);
%     gf(P+2,1:N+1,1:Q,2)   = f(P+1,:,:);
%     gf(:,:,:,2) = (P+1)*gf(:,:,:,2);
% 
%     gf(1:P+1,1:N+1,1,3)     = -f(:,:,1);
%     gf(1:P+1,1:N+1,2:Q,3)   = f(:,:,1:Q-1) - f(:,:,2:Q);
%     gf(1:P+1,1:N+1,Q+1,3)   = f(:,:,Q);
%     gf(:,:,:,3) = Q*gf(:,:,:,3);
%   %  gf = -cat(3,(N+1)*[cat(2,-D(:,1),D(:,1:N) - D(:,2:N+1),D(:,N+1)); zeros(1,N+2)],(Q+1)*[cat(1,-D(1,:),D(1:Q,:) - D(2:Q+1,:),D(Q+1,:)),zeros(Q+2,1)]);
%     
%     pU = U; 
%     pU(:,2:N+1,:,1) = pU(:,2:N+1,:,1) - gf(:,2:N+1,:,1);
%     pU(2:P+1,:,:,2) = pU(2:P+1,:,:,2) - gf(2:P+1,:,:,2);
%     pU(:,:,2:Q-1,3) = pU(:,:,2:Q-1,3) - gf(:,:,2:Q-1,3);
end

