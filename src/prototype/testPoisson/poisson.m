%Fonction pour résoudre le pbm de Poisson 
% -Df = d dans le carré [0,1]x[0,1]
%avec les conditions
%	- bL sur le bord gauche
%	- bR sur le bord droit 
%	- bU en haut 
%	- bD en bas 
%taille des conditions : 
%N+1 en haut et bas 
%Q+1 à gauche et droite
% La precision de la solution est demandée

function f = poisson(d,bL,bR,bU,bD,N,Q,tol) 
	dx = 1/N; dt = 1/Q;
	S = 0.5*(1/dx^2 + 1/dt^2)^-1;

	u = zeros(Q+1,N+1);
	u(1,:) = bU; u(end,:) = bD;
    u(:,1) = bL; u(:,end) = bR;

	ukp1 = u; 
	err = 1;
    while (err > tol) 
		for i = 2:Q
			for j = 2:N
				ukp1(i,j) = S*(d(i,j) + (u(i,j+1) + u(i,j-1))/dt^2 + (u(i+1,j) + u(i-1,j))/dx^2);
			end
		end
			err = norm(ukp1-u,2);
			u = ukp1;
    end

	f = u;
end