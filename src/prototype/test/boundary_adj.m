function U = boundary_adj(b)
    globals;
%     Bm = zeros(2*(Q+1),(N+2)*(Q+1));
%     Bm(1:Q+1,1:Q+1) = eye(Q+1);
%     Bm(Q+2:end,end-Q:end) = eye(Q+1);
%     
%     Bf = zeros(2*(N+1),(N+1)*(Q+1));
%     dia = [1 zeros(1,Q+1);zeros(1,Q+1) 1]; 
%     Bf = [];
%     for i = 1:N+1
%         Bf = blkdiag(Bf,dia);
%     end
%     
%     B = blkdiag(Bm,Bf);
%    Badj = B';
    
    U = B'*b;
end

