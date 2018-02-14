function f = gauss(mu,sigma)
    globals;
    for i = 1:N+1
        f(i) = exp(-0.5*((((i-1)/(1.0*N))-mu)/sigma).^2);
    end 
end