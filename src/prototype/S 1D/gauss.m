% Fonction gauss
% Out : 
%     - f     : une gaussienne sur [0,1] évaluée sur N + 1 points
% In : 
%     - mu    : la moyenne
%     - sigma : l'écart type 
%     - N     : le nombre d'intervals de discrétisation
% Timothéee Schmoderer 
% INSA Rouen Normandie 2017/2018

function f = gauss(mu,sigma,N)
    X = (0:N)/N;
    f = exp(-0.5*((X-mu)/sigma).^2);
end