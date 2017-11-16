% Renvoi une gaussienne normalisée dans [0,1]
% sur N+1 points
% moyenne : mu 
% varia,ce : sigma 

function f = gauss(mu,sigma,N)
    X = linspace(0,1,N+1);
    f = exp(-0.5*((X-mu)/sigma).^2);
    f = f/sum(f);
end