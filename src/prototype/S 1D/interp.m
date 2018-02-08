function V = interp(U)
    globals;
    V = zeros(Q+1,N+1,2);
    V(:,:,1) = U(1:Q+1,1:N+1,1) + U(1:Q+1,2:N+2,1);
    V(:,:,2) = U(1:Q+1,1:N+1,2) + U(2:Q+2,1:N+1,2);
    V = 0.5*V;
end

