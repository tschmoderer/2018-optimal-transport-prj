% Fonction CG 
% Out : 
%     - x   : la solution appprochée du système linéaire
%     - err : les erreurs comises par l'algorithmes à chaque itération
% In : 
%     - A : une application linéaire définie positive, donnée sous forme matricielle ou de fonctionnelle
%     - b : le second membre
% Timothée Schmoderer
% INSA Rouen Normandie 2017/2018

function x = cg(A,b)
    x = zeros(size(b));
    r = b - feval(A,x);
    p = r;
    rsold = r' * r;

    for i = 1:length(b)
        Ap = feval(A,p);
        alpha = rsold / (p' * Ap);
        x = x + alpha * p;
        r = r - alpha * Ap;
        rsnew = r' * r;
        if sqrt(rsnew) < 1e-10
              break;
        end
        p = r + (rsnew / rsold) * p;
        rsold = rsnew;
    end
end