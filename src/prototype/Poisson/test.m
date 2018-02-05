clear all
close all
clc
%% test FFT
N = 21;
Q = 19;

X       = (0:N)/N; T = (0:Q)/Q;
[XX,YY] = meshgrid(X,T); YY = flipud(YY);

normalise = @(f) f/sum(f(:)); epsilon = 1e-10;
f0 = normalise(epsilon + gauss(0.5,0.05,N));
f1 = normalise(epsilon + gauss(0.5,0.05,N));

dx  = @(m) N*(m(:,[2:end end]) - m); % dérivation selon ---> x
dt  = @(f) Q*(f([2:end end],:) - f); % dérivation selon ---> t 
dxS = @(dm) N*[-dm(:,1) , dm(:,1:end-2) - dm(:,2:end-1) , dm(:,end-1)];
dtS = @(df) Q*[-df(1,:) ; df(1:end-2,:) - df(2:end-1,:) ; df(end-1,:)];

A  = @(w)  [-dxS(w(:,:,1)) + dt(w(:,:,2)); w(end,:,2); w(1,:,2)];
AS = @(Aw) cat(3,-dx(Aw(1:Q+1,:)),dtS(Aw(1:Q+1,:)) + [Aw(Q+3,:);zeros(Q-1,N+1);Aw(Q+2,:)]);

y = [zeros(Q+1,N+1);f0;f1];

w = rand(Q+1,N+1,2);

b =  - (-dxS(w(:,:,1)) + dt(w(:,:,2)));

u = zeros(Q+1,N+1); 
u(1,:) = f1 - w(1,:,2); u(end,:) = f0 - w(end,:,2);
u(:,1) = - w(:,1,1); u(:,end) = - w(:,end,1);

dz = 1/N; dt = 1/Q;
S = 0.5*(1/dz^2 + 1/dt^2)^-1;
ukp1 = u; 
err = 1; tol = 1e-5;
k = 0;
while (err > tol) 
%    ukp1(2:end-1,2:end-1) = S*(d(2:end-1,2:end-1)+ (u(2:end-1,3:end) + u(2:end-1,1:end-2))/dt^2 + (u(3:end,2:end-1) + u(1:end-2,2:end-1))/dx^2);
 		for i = 2:Q
 			for j = 2:N
 				ukp1(i,j) = S*(b(i,j) + (u(i,j+1) + u(i,j-1))/dt^2 + (u(i+1,j) + u(i-1,j))/dz^2);
 			end
 		end
    err = norm(ukp1-u,inf);
    u = ukp1;
    k = k +1;
end

return 

pC = w + cat(3,-dx(u),dtS(u))

err = @(w) norm(A(w)-y)/norm(y);

fprintf('Error before projection: %.2e\n', err(w));
fprintf('Error after projection: %.2e\n', err(pC));


