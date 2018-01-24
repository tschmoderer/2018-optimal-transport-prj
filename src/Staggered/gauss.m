% Renvoi une gaussienne dans [0,1]x[0,1]
% sur N+1 points en X
% sur P+1 points en Y
% moyenne : muX dans la direction X, muY dans la direction Y
% variance : sigma dans la direction X et Y

function f = gauss(muX,muY,sigma,N,P)
    X = (0:N)/N;
    Y = (0:P)/P;
    [XX , YY] = meshgrid(X,Y);
    f = exp(-0.5*((XX - muX).^2 + (YY - muY).^2)/sigma^2);
end