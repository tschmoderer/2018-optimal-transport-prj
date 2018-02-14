function D = div(U)
    globals;
    D = N*(U(1:Q+1,2:N+2,1) - U(1:Q+1,1:N+1,1)) +Q*(U(2:Q+2,1:N+1,2) - U(1:Q+1,1:N+1,2));
end  