function U = interpAdj(V)
    globals; 
    U = zeros(Q+2,N+2,2);

    U(1:Q+1,1,1) = V(:,1,1);
    U(1:Q+1,2:N+1,1) = V(:,2:N+1,1) + V(:,1:N,1);
    U(1:Q+1,N+2,1) = V(:,N+1,1);

    U(1,1:N+1,2) = V(1,:,2);
    U(2:Q+1,1:N+1,2) = V(2:Q+1,:,2) + V(1:Q,:,2);
    U(Q+2,1:N+1,2) = V(Q+1,:,2);
    U = 0.5*U;
end

