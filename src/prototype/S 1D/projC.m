function pU = projC(U)
    globals;
    
    y = zeros(Q+3,N+3); y(Q+2,1:N+1) = f0; y(Q+3,1:N+1) = f1;
    
    div    = @(U) N*(U(1:Q+1,2:N+2,1) - U(1:Q+1,1:N+1,1)) +Q*(U(2:Q+2,1:N+1,2) - U(1:Q+1,1:N+1,2));
    divAdj = @(D) cat(3,N*[cat(2,-D(:,1),D(:,1:N) - D(:,2:N+1),D(:,N+1)); zeros(1,N+2)],Q*[cat(1,-D(1,:),D(1:Q,:) - D(2:Q+1,:),D(Q+1,:)),zeros(Q+2,1)]);
    A      = @(U) [[div(U) , U(1:Q+1,1,1) , U(1:Q+1,N+2,1)]; [U(Q+2,1:N+1,2) 0 0] ; [U(1,1:N+1,2) 0 0]];
    AS     = @(R) divAdj(R(1:Q+1,1:N+1)) + cat(3,[[R(1:Q+1,N+2) zeros(Q+1,N) R(1:Q+1,N+3)]; zeros(1,N+2)],[[R(Q+3,1:N+1) ; zeros(Q,N+1) ; R(Q+2,1:N+1)], zeros(Q+2,1) ]);
    bS     = @(mL,mR,fU,fD) cat(3,[mL zeros(Q+1,N) mR ; zeros(1,N+2)],[[fU ; zeros(Q,N+1) ; fD],zeros(Q+2,1)]);
    
    mL0 = zeros(Q+1,1); mR0 = zeros(Q+1,1);
    fU0 = f1; fD0 = f0;
    
    B = bS(mL0 - U(1:Q+1,1,1),mR0 - U(1:Q+1,N+2,1),fU0 - U(1,1:N+1,2),fD0 - U(Q+2,1:N+1,2));
    
    w  = div(U + B);
    x1 = poisson2d_Neumann(-w);
    
    pU = U - divAdj(x1) + B;
end

%     x = zeros(Q+3,N+3);
%     b = y - A(U);
%     r = b - A(AS(x));
%     p = r;
%     rold = sum(r(:).*r(:));
% 
%     for i = 1:(Q+3)*(N+3)
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
%     pU = U + AS(x);

%    U(1:Q+1,[1 end],1) = 0; U([1 end],1:N+1,2) = [f1 ; f0]; 


% function D = div(U)
%     globals; 
% 	D = zeros(Q+1,N+1);
% 	D = N*(U(1:Q+1,2:N+2,1) - U(1:Q+1,1:N+1,1)) +Q*(U(2:Q+2,1:N+1,2) - U(1:Q+1,1:N+1,2));
% end

% function U = divAdj(D)
%     globals;
%     U = zeros(Q+2,N+2,2);
%     
%     U(1:Q+1,1,1) = -D(:,1);
%     U(1:Q+1,2:N+1,1) = D(:,1:N) - D(:,2:N+1);
%     U(1:Q+1,N+2,1) = D(:,N+1);
%     U(:,:,1) = N*U(:,:,1);
% 
%     U(1,1:N+1,2) = -D(1,:);
%     U(2:Q+1,1:N+1,2) = D(1:Q,:) - D(2:Q+1,:);
%     U(Q+2,1:N+1,2) = D(Q+1,:);
%     U(:,:,2) = Q*U(:,:,2);
% end

% function Au = A(U)
%     globals;
%     Au = zeros(Q+3,N+3);
%     
%     Au(1:Q+1,1:N+1) = div(U);
%     Au(1:Q+1,N+2)   = U(1:Q+1,1,1);   % frontieres de mbar
%     Au(1:Q+1,N+3)   = U(1:Q+1,N+2,1); 
%     Au(Q+2,1:N+1)   = U(Q+2,1:N+1,2);   % frontieres de fbar
%     Au(Q+3,1:N+1)   = U(1,1:N+1,2);
% end

% function U = AS(R)
%     globals;
%     U = zeros(Q+2,N+2,2);
%     
%     U(1:Q+2,1:N+2,:) = divAdj(R(1:Q+1,1:N+1));
%     U(1:Q+1,1,1) = U(1:Q+1,1,1) + R(1:Q+1,N+2);
%     U(1:Q+1,N+2,1) = U(1:Q+1,N+2,1) + R(1:Q+1,N+3);
%     U(1,1:N+1,2) = U(1,1:N+1,2) + R(Q+3,1:N+1);
%     U(Q+2,1:N+1,2) = U(Q+2,1:N+1,2) + R(Q+2,1:N+1);
% end