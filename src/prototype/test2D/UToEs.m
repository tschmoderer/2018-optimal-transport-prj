function [mbar,fbar] = UToEs(U,N,P,Q)
    mbar = reshape(U(1:2*(N+2)*(P+2)*(Q+1)),[P+2,N+2,Q+1,2]);
    fbar = reshape(U(2*(N+2)*(P+2)*(Q+1)+1:end),[P+1,N+1,Q+2]);
end

