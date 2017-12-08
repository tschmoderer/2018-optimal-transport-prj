% Use of FFT for FD approximation to u for
%
%      u_xx + u_yy  = 4,  0 < x,y < 1
%      u(0,y) = y^2, u(1,y) = 1 + y^2, u(x,0) = x^2, u(x,1) = 1 + x^2
%
%  with h = 1/8. This eqn has exact solution uex = x^2 + y^2
%  For simplicity, this script uses FFTs in both directions, rather than 
%  an FFT in one direction and a tridiagonal solver in the other.
%  Note that u_{i,j} corresponds to ROW i and COLUMN j of a matrix of u's;
%  to see the gridpoint values in their correct geometrical orientation
%  so increasing x runs to the right and increasing y upward, 
%  you must display u(:,m:-1:1)'.  Note that in this case, the FD method
%  has zero truncation error since u is parabolic, so the FD soln is exact.
%  Script makes a surface plot of u at the interior grid points.

function u = poissonFFT(d,bL,bR,bU,bD,N,P)
  % m =  Number of interior gridpoints in original domain
  L = 1;        % Domain size
  Ne = 2*N+2;   % Extended-domain gridpoints...should be a power even
  NeX =  2*N+2; NeT = 2*P+2;
  h = 2*L/Ne;   % Grid spacing
  hX = 2*L/NeX; hT = 2*L/NeT;
  x = h*(1:N);
  y = h*(1:N);
  x= hX*[1:N];
  y = hT*[1:P];
  
  f = -d;
  
  % Correct for Dirichlet BCs
  u0y = bD; % y = 0
  u1y = bU; % y = 1
  ux0 = bL; % x = 0
  ux1 = bR; % x = 1

  f(:,1) = f(:,1) - ux0'/hT^2;
  f(:,end) = f(:,end) - ux1'/hT^2;
  f(1,:) = f(1,:) - u0y/hX^2;
  f(end,:) = f(end,:) - u1y/hX^2;
  
  g = [zeros(N,1) f zeros(N,1) -f(:,N:-1:1)];    %  Odd extension of f in x
  g = [zeros(1,Ne); g; zeros(1,Ne); -g(N:-1:1,:)]; %  Odd extension of g in y
  ghat = fft2(g);
  [I J]  = meshgrid(0:(Ne-1),0:(Ne-1));
  mu = (4/h^2)*(sin(I*pi/Ne).^2 + sin(J*pi/Ne).^2);
  mu(1,1) = 1;   % Avoid 0/0 division; vhat(1,1) is known a priori to be 0
  v = real(ifft2(-ghat./mu));
  u = v(1:(N+2),1:(N+2)); % Extract u(i,j) from its odd extension
  
%  Plot out solution in interior and print out max-norm error

%   [X, Y] = meshgrid(x,y);
%   uex = X.^2 + Y.^2;
%   err = abs(u(2:(m+1),2:(m+1)) - uex); 
%   err = norm(err(:),inf)
%   surf(X,Y,u(2:(m+1),2:(m+1)))
%   xlabel('x')
%   ylabel('y')
%   zlabel('u')
%   title('FD/FFT for 2D Poisson Eqn')
 u = u(2:(N+1),2:(N+1));
end
