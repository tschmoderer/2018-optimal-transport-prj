% renvoi l'indicatrice de [a,b]x[c,d]
% ou 0<=a<b<=1 && 0<=c<d<=1
% avec [0,1]x[0,1] discrétisé dans (N+1)(P+1) pts

function f = indicatrix(a,b,c,d,N,P)    
    X = (0:N)/N;
    Y = (0:P)/P;
    [XX , YY] = meshgrid(X,Y);

    f = zeros(size(XX));
    f(a<=XX & XX<=b & c<=YY & YY<=d) = 1;
end

