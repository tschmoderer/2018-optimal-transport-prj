
function d = divergence(U)
    globals;
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
    
    d = D*U;
    
end

