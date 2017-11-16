% Renvoi une gaussienne normalisée dans [0,1]x[0,1]
% sur (N+1)x(P+1) points
% moyennes : muX,muY
% variances : sigmaX, sigmaY 
% plancher minimal : min

function f = gauss(muX,muY,sigmaX,sigmaY,N,P,min)
    [X, Y] = meshgrid(linspace(0,1,N+1),linspace(0,1,P+1));
    f = min + exp(-0.5*((X-muX).^2/sigmaX^2 + (Y-muY).^2/sigmaY^2));
    f = f/sum(f(:));
end