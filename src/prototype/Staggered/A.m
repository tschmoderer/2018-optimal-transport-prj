function Au = A(U)
    globals;
    Au = zeros(P+3,N+3,Q+3);
    Au(1:P+1,1:N+1,1:Q+1) = div(U);

    Au(1:P+1,N+2,1:Q+1)   = U(1:P+1,N+2,1:Q+1,1);   % frontieres de mbar1
    Au(1:P+1,N+3,1:Q+1)   = U(1:P+1,1  ,1:Q+1,1); 

    Au(P+2,1:N+1,1:Q+1)   = U(P+2,1:N+1,1:Q+1,2);   % frontieres de mbar 2
    Au(P+3,1:N+1,1:Q+1)   = U(1  ,1:N+1,1:Q+1,2);
    
    Au(1:P+1,1:N+1,Q+2)   = U(1:P+1,1:N+1,Q+2,3);   % frontieres de fbar
    Au(1:P+1,1:N+1,Q+3)   = U(1:P+1,1:N+1,1  ,3);
end

