function R = cost(V)
    globals;
    R = (V(1:(N+1)*(P+1)*(Q+1)).^2 + V((N+1)*(P+1)*(Q+1)+1:2*(N+1)*(P+1)*(Q+1)).^2);
    R = R./V(2*(N+1)*(P+1)*(Q+1):end);
    R = sum(R);
end

