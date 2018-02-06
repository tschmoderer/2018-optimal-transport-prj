% Solve -Delta(u) = f
% on ]0,1[x]0,1[
% N+1 pts en x, Q+1 pts en t 
% with u'(x=0,:) = u'(1,:) = 0
% with u(:,t=0) = f1, u(:,1) = f0

function u = poisson_neumann_2d(f,f0,f1,N,Q)
	S = 0.5*(N^2 + Q^2)^-1; R = 1.;
    u = zeros(Q+1,N+3); u1 = u;
    u(1,2:N+2) = f1; u(Q+1,2:N+2) = f0;
    % f is (Q+1)x(N+1)
    while R > 1e-5
       for i = 2:Q
           for j = 2:N+2
            u1(i,j) = S*(f(i,j-1) + (u(i,j+1) + u(i,j-1))*Q^2 + (u(i+1,j) + u(i-1,j))*N^2);
%           u1(i,j) = 0.25*(h^2*f(i,j) + u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1));
           end
       end
      
        u1(1,2:N+2) = f1; u1(Q+1,2:N+2) = f0;
        u1(:,1) = u(:,2); u1(:,N+3) = u(:,N+2);
        u1(:,1) = 0; u1(:,N+3) = 0;
        
        R = norm(u(:,2:N+2)-u1(:,2:N+2))/norm(u(:,2:N+2))
        u(:,2:N+2) = u1(:,2:N+2);
    end
    u = u(:,2:N+2);
end

