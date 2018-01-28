% Fonction CG 
% Out : 
%     - x   : la solution appprochée du système linéaire
%     - err : les erreurs comises par l'algorithmes à chaque itération
% In : 
%     - A : une application linéaire définie positive, donnée sous forme matricielle ou de fonctionnelle
%     - b : le second membre
% Timothée Schmoderer
% INSA Rouen Normandie 2017/2018

function [x, err] = cg(A,b)
    dotp = @(x,y) sum(x(:).*y(:));
    
    % initial guess
    x0 = zeros(size(b));
    
    if isnumeric(A)
        r0 = b - A*x0;
    else
        r0 = b - feval(A,x0);
    end
    p0 = r0;
    
    k = 1;
    niter = 2000;
    err = zeros(niter,1);
        
    while k < niter
       if isnumeric(A)
           alpha = dotp(r0,r0)/dotp(p0,A*p0);
           r1    = r0 - alpha*A*p0;
       else
           alpha =  dotp(r0,r0)/dotp(p0,feval(A,p0));
           r1    =  r0 - alpha*feval(A,p0);
       end
       x1    = x0 + alpha*p0;
       
       err(k) = norm(r1,1);
       if err(k) < 1e-9
           break
       end
       
       beta = dotp(r1,r1)/dotp(r0,r0);
       p1   = r1 + beta*p0; 
       
       r0 = r1;
       p0 = p1;
       x0 = x1;
       k = k+1;
    end
    if k == niter 
        warning('CG converged without reaching the prescribed precision');
    end
    err = err(1:k);
    x = x1;
end

