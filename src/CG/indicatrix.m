% Fonction indicadrix
% Out : 
%     - f = 1 si x \in [a,b], 0 sinon
% In  :
%     - a : borne à gauche 
%     - b : borne à droite 
%     - N : nombre d'intervals de subdivision
% /!\ 0<= a < b <= 1
% Timothée Schmoderer
% INSA Rouen Normandie 2017/2018

function f = indicatrix(a,b,N)
    x = (0:N)/N;
    f = zeros(size(x));
    f(a<=x & x<=b) = 1;
end

