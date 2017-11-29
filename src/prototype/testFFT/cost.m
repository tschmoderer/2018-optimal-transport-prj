function R = cost(V)
    globals;
    m = reshape(V(1:(N+1)*(Q+1)),Q+1,N+1);
    f = reshape(V((N+1)*(Q+1)+1:end),Q+1,N+1);
    
    R = 0.5*m.^2./f;
    
   % idx = find(f == 0 & m == 0);
   % R(idx) = 0;
   % idx = find(f < 0);
   % R(idx) = 10^5;
    R = sum(R(:));
%     R = 0.5*V(1:(N+1)*(Q+1)).^2./V((N+1)*(Q+1)+1:end);
%     R = sum(R);
end

