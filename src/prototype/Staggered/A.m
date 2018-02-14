function Au = A(U)
    globals;
    Au = zeros(Q+3,N+3);
    Au(1:Q+1,1:N+1) = div(U);

    Au(1:Q+1,N+2)   = U(1:Q+1,1,1);   % frontieres de mbar
    Au(1:Q+1,N+3)   = U(1:Q+1,N+2,1); 

    Au(Q+2,1:N+1)   = U(Q+2,1:N+1,2);   % frontieres de fbar
    Au(Q+3,1:N+1)   = U(1,1:N+1,2);
end

