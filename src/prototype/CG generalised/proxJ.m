% Validée
% generealised cost 

function Pw = proxJ(w,b,g)
    mt = w(:,:,1);
    ft = w(:,:,2);
    
    x0 = 1000;
    x1 = 2000;
    k = 0;

    % Newton
    while norm(x0-x1,1) > 1e-5 && k < 1500
        x0 = x1;
        poly = x0.^(1-b).*(x0-ft).*((x0.^b+g).^2)-0.5*g*b*mt.^2;
        dpoly = (1-b)*x0.^(-b).*(x0-ft).*((x0.^b+g).^2) + x0.^(1-b).*((x0.^b+g).^2 +2*b*(x0-ft).*x0.^(b-1).* (x0.^b+g) );
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


