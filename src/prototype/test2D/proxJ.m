function pJ = proxJ(V,g)
    globals;
    
    mt1 = V(1:(N+1)*(P+1)*(Q+1));
    mt2 = V((N+1)*(P+1)*(Q+1)+1:2*(N+1)*(P+1)*(Q+1));
    ft = V(2*(N+1)*(P+1)*(Q+1):end);
    
    
    mt1 = reshape(mt1,[P+1,N+1,Q+1]);
    mt2 = reshape(mt2,[P+1,N+1,Q+1]);
    ft = reshape(ft,[P+1,N+1,Q+1]);
      
    Poly = @(x)(-g*(mt1.^2 + mt2.^2) + (x-ft).*(x+2*g).^2);
	dP = @(x)(2*(x+2*g).*(x-ft)+(x+2*g).^2);
 
	x0 = 1000;
	x1 = 0;
	k = 0;
	% Newton
	while norm(x0-x1,1) > 1e-5 && k < 50
		x0 = x1;
		x1 = x0 - Poly(x0)./dP(x0);
		k = k+1;
	end

	Pf = x1;
	Pm = x1.*mt./(x1 + 2*g);
	idx = find(x1 <= 0);
	Pm(idx) = 0;
	Pf(idx) = 0;

    pJ = [reshape(Pm,(N+1)*(Q+1),1);reshape(Pf,(N+1)*(Q+1),1)];
end

