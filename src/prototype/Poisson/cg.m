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
    r = b - A(x);
    p = r;
    rsold = sum(r(:).*r(:));

    for i = 1:prod(size(b))
        Ap = A(p);
        alpha = rsold / sum(p(:).*Ap(:));
        x = x + alpha * p;
        r = r - alpha * Ap;
        rsnew = sum(r(:).*r(:));
        if sqrt(rsnew) < 1e-10
              break;
        end
        p = r + (rsnew / rsold) * p;
        rsold = rsnew;
    end
end