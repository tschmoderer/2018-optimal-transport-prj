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
	u(1,:) = bU; u(end,:) = bD; %!!%
    %u(:,1) = bL; u(:,end) = bR;

%     % !! 
%     u(1,:) = d(1,:) - bU; u(end,:) = d(end,:) - bD;
%     u(:,1) = d(:,1) - bL; u(:,end) = d(:,end) - bR;  
%     % !! 
%     
	ukp1 = u; 
	err = 1;
    while (err > tol) 
        ukp1(2:end-1,2:end-1) = S*(d(2:end-1,2:end-1)+ (u(2:end-1,3:end) + u(2:end-1,1:end-2))/dt^2 + (u(3:end,2:end-1) + u(1:end-2,2:end-1))/dx^2);
% 		for i = 2:Q
% 			for j = 2:N
% 				ukp1(i,j) = S*(d(i,j) + (u(i,j+1) + u(i,j-1))/dt^2 + (u(i+1,j) + u(i-1,j))/dx^2);
% 			end
% 		end
		err = norm(ukp1-u,inf);
		u = ukp1;
    end

	f = u;















% % Solving the 2-D Poisson equation by the Finite Difference
% ...Method 
% % Numerical scheme used is a second order central difference in space
% ...(5-point difference)
% 
% %%
% %Specifying parameters
%  
% nx = N+1; ny = Q+1;
% niter=1000;                      %Number of iterations 
% dx = 1/N;                     %Width of space step(x)
% dy = 1/Q;                     %Width of space step(y)
% x = 0:dx:1;                        %Range of x(0,2) and specifying the grid points
% y = 0:dy:1;                        %Range of y(0,2) and specifying the grid points
% b = zeros(nx,ny);                  %Preallocating b
% pn = zeros(nx,ny);                 %Preallocating pn
% 
% %%
% % Initial Conditions
% p=zeros(nx,ny);                  %Preallocating p
% 
% %%
% %Boundary conditions
% p(:,1) = bL;
% p(:,ny) = bR;
% p(1,:) = bU;                  
% p(nx,:) = bD;
% 
% %%
% %Source term
% % b(round(ny/4),round(nx/4))=3000;
% % b(round(ny*3/4),round(nx*3/4))=-3000;
% b = d;
% %%
% i=2:nx-1;
% j=2:ny-1;
% %Explicit iterative scheme with C.D in space (5-point difference)
% for it=1:niter
%     pn=p;
%     p(i,j)=((dy^2*(pn(i+1,j)+pn(i-1,j)))+(dx^2*(pn(i,j+1)+pn(i,j-1)))-(b(i,j)*dx^2*dy*2))/(2*(dx^2+dy^2));
%     %Boundary conditions 
% %     p(:,1)=0;
% %     p(:,ny)=0;
% %     p(1,:)=0;                  
% %     p(nx,:)=0;
% end
% 
% %%
% %Plotting the solution
% h=surf(x,y,p,'EdgeColor','none');       
% shading interp
% %axis([-0.5 2.5 -0.5 2.5 -100 100])
% title({'2-D Poisson equation';['{\itNumber of iterations} = ',num2str(it)]})
% xlabel('Spatial co-ordinate (x) \rightarrow')
% ylabel('{\leftarrow} Spatial co-ordinate (y)')
% zlabel('Solution profile (P) \rightarrow')
% 
% f = p;







end