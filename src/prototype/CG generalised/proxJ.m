% Validée
% generealised cost 
% with obstacle is matrix [Q+1]x[N+1] 
% obstacle = 1 > 0 si ok 
%          = +\infty si transport impossible

function Pw = proxJ(w,b,g,obstacle)
    mt = w(:,:,1);
    ft = w(:,:,2);
    
    x0 = 1000;
    x1 = 2000;
    k = 0;

    g = g*obstacle;
    % Newton
    while norm(x0-x1,1) > 1e-5 && k < 1500
        x0 = x1;
        if b == 1 % cas transport
            poly = (x0-ft).*((x0+g).^2)-0.5*g.*mt.^2;
            dpoly = 2*(x0+g).*(x0-ft)+(x0+g).^2;
        elseif b == 0 % cas interpolation 
            x0    = ft;
            poly  = 0; 
            dpoly = 1;
        else
            poly = x0.^(1-b).*(x0-ft).*((x0.^b+g).^2)-0.5*b*g.*mt.^2;
            dpoly = (1-b)*x0.^(-b).*(x0-ft).*((x0.^b+g).^2) + x0.^(1-b).*((x0.^b+g).^2 +2*b*(x0-ft).*x0.^(b-1).* (x0.^b+g) );
        end
        x1 = x0 - poly./dpoly;
        k = k+1;
    end
    Pf = x1;
    Pm = (Pf.^b).*mt./(Pf.^b + g);

    idx = find(Pf <= 0);
    Pm(idx) = 0;
    Pf(idx) = 0;
    
    Pw = zeros(size(w));
    Pw(:,:,1) = real(Pm);
    Pw(:,:,2) = real(Pf);
end


