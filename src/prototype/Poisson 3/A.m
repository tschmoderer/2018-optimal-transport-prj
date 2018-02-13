function Au = A(U)
    globals;
    Au = zeros(P+3,N+3,Q+2);
    Au(1:P+1,1:N+1,1:Q) = div(U);

    Au(1:P+1,N+2,1:Q)   = U(1:P+1,1,1:Q,1); % frontiere de mbar1
    Au(1:P+1,N+3,1:Q)   = U(1:P+1,N+2,1:Q,1);
		
    Au(P+2,1:N+1,1:Q)   = U(1,1:N+1,1:Q,2); % frontiere de mbar2
    Au(P+3,1:N+1,1:Q)   = U(P+2,1:N+1,1:Q,2);
		
    Au(1:P+1,1:N+1,Q+1) = U(1:P+1,1:N+1,1,3); % frontieres de fbar
    Au(1:P+1,1:N+1,Q+2) = U(1:P+1,1:N+1,Q+1,3);
end

