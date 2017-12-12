function [Pm,Pf] = proxJ(mt,ft)
	x0 = 1000;
	x1 = 0;
	k = 0;
    g = 1.0;
	% Newton
	while norm(x0-x1,1) > 1e-5 && k < 15000
		x0 = x1;
        poly = (x0-ft).*((x0+g).^2)-0.5*g*mt.^2;
        dpoly = 2*(x0+g).*(x0-ft)+(x0+g).^2;
        ddpoly = 2*(3*x0+2*g-ft);
        x1 = x0 - 2*poly.*dpoly./(2*dpoly.^2 - poly.*ddpoly);
        k = k+1;
    end
   k
	Pf = x1;
	Pm = Pf.*mt./(Pf + g);
    
 	idx = find(Pf <= 0);
 	Pm(idx) = 0;
 	Pf(idx) = 0;
end


