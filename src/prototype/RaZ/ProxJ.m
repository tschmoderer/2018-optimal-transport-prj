function [Pm,Pf] = ProxJ(m,f,gamma)
    P = @(X)((X-f).*(X+2*gamma).^2 - gamma*m.^2);
    dP = @(X)(2*(X-f).*(X+2*gamma)+(X+2*gamma).^2);
    
    x1 = 1000; x0 = 0; k = 0;
    
    while norm(x1-x0,1) > 1e-5 && k < 50
        x0 = x1;
        x1 = x0 - P(x0)./dP(x0);
        k = k+1;
    end
    
    Pm = x1.*m./(x1+2*gamma);
    Pf = x1;
    
    idx = find(x1 <= 0);
    Pm(idx) = 0;
    Pf(idx) = 0;
end