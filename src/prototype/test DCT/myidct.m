function df = myidct(f,N)
    df= zeros(size(f));
    a = sqrt(2/(N+1));
    b = ones(1,N+1); b(1) = 1/sqrt(2);
    for k = 1:N+1
%         s = 0;
%         for n = 1:N+1
%            s = s + f(n)*b(n)*cos(pi*(n-1)*(2*k-1)/(2*(N+1))); 
%         end
        s = sum(f.*b.*cos(pi*(0:N)*(2*k-1)/(2*(N+1))));
        df(k) = a*s;
    end
end

