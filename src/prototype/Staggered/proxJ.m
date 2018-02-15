% Fonction proxJ 
% Out : 
%     - Pw : l'opérateur proximale associé à J évalué en w
% In  : 
%     - w  : le point d'évaluation 
%     - b  : le coefficient beta du cout généralisé
%     - g  : le coefficient gamma qui pondère l'opérateur proximal 
%     - obstacle : la matrice des obstacle
% Timothée Schmoderer 
% INSA Rouen Normandie 2017/2018

function Pw = proxJ(w,b,g,obstacle)
    globals; 

    mt = w(:,:,:,1:2); ft = w(:,:,:,3);
    x0 = ones(size(ft)); x1 = 2*ones(size(ft)); k = 0;

    while (max(sum(sum(abs(x0-x1)))) > 1e-5  && k < 1500)
        x0 = x1;
        if (b == 1)  
            poly  = (x0-ft).*(x0+g).^2 - 0.5*g*(mt(:,:,:,1).^2 + mt(:,:,:,2).^2);
            dpoly = 2*(x0+g).*(x0-ft) + (x0+g).^2;
        elseif (b == 0)
            x1 = ft;
            break
        else 
            poly = x0.^(1.0-b).*(x0-ft).*((x0.^b+g).^2)-0.5*b*g*(mt(:,:,:,1).^2 + mt(:,:,:,2).^2);
            dpoly = (1.0-b)*x0.^(-b).*(x0-ft).*((x0.^b+g).^2) + x0.^(1-b)*((x0.^b+g).^2 +2*b*(x0-ft).*x0.^(b-1).*(x0.^b+g) );
        end
        
        idx = find(x0 > eps); x1(idx) = x0(idx) - poly(idx)./dpoly(idx);
        idx = find(x0 < eps); x1(idx) = eps;
    
        k = k+1;
    end 
    idx = find(x1 < eps | obstacle > 0); x1(idx) = eps;
    Pw = zeros(size(w));
    Pw(:,:,:,3) = x1;
    Pw(:,:,:,1:2) = (x1.^b).*mt./(x1.^b+g);
end


