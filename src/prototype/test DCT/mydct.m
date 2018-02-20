% f is N+1
function df = mydct(f,N)   
    df= zeros(size(f));
    a = sqrt(2/(N+1));
    b = ones(1,N+1); b(1) = 1/sqrt(2);
    for k = 1:N+1
        s = sum(f.*cos(pi*((1:N+1) - 0.5)*(k-1)/(N+1)))*b(k);   
        df(k) = a*s;
    end   
end


