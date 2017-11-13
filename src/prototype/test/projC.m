function pC = projC(U)
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
%   Badj = B';
    
%     Dm = zeros((N+1)*(Q+1),(N+2)*(Q+1));
%     for i = 1:(N+1)*(Q+1)
%         for j = 1:(N+2)*(Q+1)
%             if i == j 
%                 Dm(i,j) = -1;
%             elseif j == i+Q+1
%                 Dm(i,j) = 1;
%             end
%         end
%     end
%     dia = zeros(Q+1,Q+2);
%     for i = 1:Q+1
%         for j = 1:Q+2
%             if i == j 
%                 dia(i,j) = -1;
%             elseif j == i+1
%                 dia(i,j) = 1;
%             end
%         end
%     end
%     Df = [];
%     for i = 1:N+1
%         Df = blkdiag(Df,dia);
%     end
%     
%     D = [Dm Df];
%   Dadj = D';
    

    
    
    pC = P*U + Cst;
end

