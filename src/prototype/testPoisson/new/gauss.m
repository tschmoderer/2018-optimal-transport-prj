% Renvoi une gaussienne dans [0,1]
% sur N+1 points
% moyenne : mu 
% variance : sigma 

function f = gauss(mu,sigma,N)
    X = [0:N]/N;
    f = exp(-0.5*((X-mu)/sigma).^2);
end