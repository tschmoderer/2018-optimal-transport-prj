% Renvoi une gaussienne dans [0,1]x[0,1]
% sur N+1 points en X
% sur P+1 points en Y
% moyenne : muX dans la direction X, muY dans la direction Y
% variance : sigmaX dans la direction X, sigmaY dans la direction Y

function f = gauss(muX,muY,sigmaX,sigmaY,N,P)
    X = [0:N]/N;
    Y = [0:P]/P;
    [XX , YY] = meshgrid(X,Y);
    f = exp(-0.5*(((XX-muX)/sigmaX).^2 + ((YY-muY)/sigmaY).^2));
end