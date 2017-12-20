% renvoi l'indicatrice de [a,b]
% ou 0<=a<b<=1
% avec [0,1]discrétisé dan N+1 pts

function f = indicatrix(a,b,N)
    x = [0:N]/N;
    f = zeros(size(x));
    f(find(a<=x & x<=b)) = 1;
end

