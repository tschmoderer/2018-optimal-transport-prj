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
    g = g*obstacle;
    
    mt = w(:,:,1);
    ft = w(:,:,2);
    
    x0 = 1;
    x1 = 2;
    k = 0;

    % Newton
    while norm(x0-x1,1) > 1e-5 && k < 1500
        x0 = x1;
        if b == 1 % cas transport
            poly = (x0-ft).*((x0+g).^2)-0.5*g.*mt.^2;
            dpoly = 2*(x0+g).*(x0-ft)+(x0+g).^2;
        elseif b == 0 % cas interpolation 
            x1 = ft;
            break
        else
            poly = x0.^(1-b).*(x0-ft).*((x0.^b+g).^2)-0.5*b*g.*mt.^2;
            dpoly = (1-b)*x0.^(-b).*(x0-ft).*((x0.^b+g).^2) + x0.^(1-b).*((x0.^b+g).^2 +2*b*(x0-ft).*x0.^(b-1).* (x0.^b+g) );
        end
        x1 = x0 - poly./dpoly;
        k = k+1;
    end
    
    Pf = x1;
    Pm = (Pf.^b).*mt./(Pf.^b + g);
% 
%     idx = find(Pf <= 0 | obstacle > 0); % ou la racine est négative ou la ou il y a obstacles
%     Pm(idx) = 0;
%     Pf(idx) = epsilon;
%     
    Pw = zeros(size(w));
    Pw(:,:,1) = real(Pm);
    Pw(:,:,2) = real(Pf);
    
%     Pw(:,:,2) = x1;
%     Pw(:,:,1) = x1.*mt./(x1+g);
    
    for i = 1:Q+1
        for j = 1:N+1
            if x1(i,j) <= 0
                Pw(i,j,1) = 0;
                Pw(i,j,1) = 1e-10;
            end
        end
    end   
end


