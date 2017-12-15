% résoud l'équation Ax = b par la méthode du gradient conjugué
% où A est une fonction

function [x] = cg(A,b)
    % initial guess
    x0 =  ones(size(b));
    r0 = b - feval(A,x0);
    p0 = r0;
    k = 0;

    while k < 1000
       alpha =  (sum(sum(sum(r0'*r0))))/(sum(sum(sum(p0'*feval(A,p0)))));
       x1    = x0 + alpha*p0;
       r1    = b - feval(A,x1);
       
       if norm(r1,1) < 1e-9
           break
       end
       
       beta = (sum(sum(sum(r1'*r1))))/(sum(sum(sum(r0'*r0))));
       p1   = r1 + beta*p0; 
       
       r0 = r1;
       p0 = p1;
       x0 = x1;
       k = k+1;
    end
    k
    x = x1;
end

