function Pw = proxJ(w,g,N,Q)
   
    x0 = 1000; 
    x1 = 2000;
    mt = w(:,:,1); 
    ft = w(:,:,2);
    k = 0;

    while norm(x0-x1,1) > 1e-10 && k < 1500
        x0 = x1;
        poly  = (x0 - ft).*(x0 + g).*(x0 + g)-0.5*g*(mt.*mt);
		dpoly = 2.0*(x0 + g).*(x0 - ft) + (x0 + g).*(x0 + g);
     %   poly  = (x0 - ft).*((x0 + g).^2)-0.5*g*(mt.^2);
      %  dpoly = 2.0*(x0 + g).*(x0 - ft) + (x0 + g).^2;
    % 			ddpoly = 2.0*(2.0*x0-ft+g) + 2.0*(x0+g)
        x1 = x0 - poly./dpoly;
    % 			x1 = x0 - 2*poly*dpoly/(2*dpoly**2-poly*ddpoly)
        k = k+1;
    end


     Pw = zeros(size(w));
     Pw(:,:,2) = x1;
     Pw(:,:,1) = x1.*mt./(x1+g);
     for i = 1:Q+1
         for j = 1:N+1
             if x1(i,j) <0 
              Pw(i,j,1) = 0;
              Pw(i,j,2) = 0;
             end
         end
     end
end

