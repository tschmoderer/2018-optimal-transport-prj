function f = gauss(muX,muY,sigma)
    globals;
    f = zeros(N+1,P+1);
    for i = 1:N+1
        for j = 1:P+1
            x = (i-1)/(1.0*N);
            y = (j-1)/(1.0*P);
            f(i,j) = exp(-0.5*(((x-muX)^2 + (y-muY)^2)/sigma^2));
        end 
    end
end