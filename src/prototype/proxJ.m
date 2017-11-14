function pJ = proxJ(V,g)
    globals;
    mt = reshape(V(1:(N+1)*(Q+1)),Q+1,N+1);
    ft = reshape(V((N+1)*(Q+1)+1:end),Q+1,N+1);
    
    Poly = @(x)(-0.5*g*mt.^2+(x-ft).*(x+g).^2);
	dP = @(x)(2*(x+g).*(x-ft)+(x+g).^2);
 
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
	Pm = x1.*mt./(x1 + g);
	idx = find(x1 <= 0);
	Pm(idx) = 0;
	Pf(idx) = 0;

    pJ = [reshape(Pm,(N+1)*(Q+1),1);reshape(Pf,(N+1)*(Q+1),1)];
end

