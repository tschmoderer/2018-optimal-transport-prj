function f = poisson_neumann(d,bL,bR,bU,bD,N,Q,tol) 
    dx = 1/N; dt = 1/Q;
    S = 1/(2*dx^2+2*dt^2);

    m2 = bR; % G2
    m4 = bL; % G4

    u = zeros(Q+1,N+1);

    u(1,:) = bU;
    u(end,:) = bD;

    ukp1 = u; 
    err = 1;
    while (err > tol) 
            for i = 2:Q
                for j = 2:N
                    ukp1(i,j) = S*(d(i,j)+dt^2*(u(i-1,j)+u(i+1,j))+dx^2*(u(i,j-1)+u(i,j+1)));
                    ukp1(i,1) = ukp1(i,2) - dx*m2(i);
                    ukp1(i,N+1) = ukp1(i,N) + dx*m4(i);
                end
            end
        err = norm(ukp1-u,inf);
        u = ukp1;
    end

    f = u;
end