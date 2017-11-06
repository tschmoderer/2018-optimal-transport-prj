function [Pm, Pf] = proxJ(mt,ft,g)
	P = @(x)(-g*mt.^2+(x-ft).*(x+2*g).^2);
	dP = @(x)(2*(x+2*g).*(x-ft)+(x+2*g).^2);
    ddP = @(x)(4*(x+2*g)+2*(x-ft));
    
	x0 = 1000;
	x1 = 0;
	k = 0;
	% Newton
	while norm(x0-x1,1) > 1e-5 && k < 50
		x0 = x1;
		x1 = x0 - P(x0)./dP(x0);
		k = k+1;
	end

	Pf = x1;
	Pm = x1.*mt./(x1 + 2*g);
	idx = find(x1 <= 0);
	Pm(idx) = 0;
	Pf(idx) = 0;
end