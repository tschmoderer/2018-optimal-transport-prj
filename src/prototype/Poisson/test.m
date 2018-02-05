clear all
close all
clc
%% test FFT
N = 5;
Q = 5;

X       = (0:N)/N; T = (0:Q)/Q;
[XX,YY] = meshgrid(X,T); YY = flipud(YY);

normalise = @(f) f/sum(f(:)); epsilon = 1e-10;
f0 = normalise(epsilon + gauss(0.5,0.05,N));
f1 = normalise(epsilon + gauss(0.5,0.05,N));

dx  = @(m) N*(m(:,[2:end end]) - m); % dérivation selon ---> x
dt  = @(f) Q*(f([2:end end],:) - f); % dérivation selon ---> t 
dxS = @(dm) N*[-dm(:,1) , dm(:,1:end-2) - dm(:,2:end-1) , dm(:,end-1)];
dtS = @(df) Q*[-df(1,:) ; df(1:end-2,:) - df(2:end-1,:) ; df(end-1,:)];

A  = @(w)  [-dxS(w(:,:,1)) + dt(w(:,:,2)); w(end,:,2) ;w(1,:,2)];
AS = @(Aw) cat(3,-dx(Aw(1:Q+1,:)),dtS(Aw(1:Q+1,:)) + [Aw(Q+3,:);zeros(Q-1,N+1);Aw(Q+2,:)]);

y = [zeros(Q+1,N+1);f0;f1];

w = rand(Q+1,N+1,2);
pC = w + AS(poisson(y - A(w),zeros(1,Q+1),zeros(1,Q+1),f1,f0,N,Q,1e-5))



% f = -dxS(w(:,:,1)) + dt(w(:,:,2));
% hx                    = 1/(N+1);
% dn                    = 0:1:N;
% depn                  = 2*cos(pi*dn/(N+1))-2;
% 
% hy                    = 1/(Q+1);
% dm                    = 0:1:Q;
% depm                  = 2*cos(pi*dm/(Q+1))-2;
% 
% denom                = repmat(depn(:)/hx^2,1,Q+1) + repmat(depm(:)'/hy^2,N+1,1);
% denom(denom(:)==0)  = 1;
% 
% 
% fhat                  = dct2(f);
% uhat                  = -(fhat)./denom';
% 
% res = idct2(uhat);

% pC = w0 + cat(3,-dx(res),dtS(res));
% 
% %% check div=0
% y = [zeros(Q+1,N+1);f0;f1];
% err = @(w) norm(A(w)-y)/norm(y);
% 
% % fprintf('Error before projection: %.2e\n', err(w));
% % fprintf('Error after projection: %.2e\n', err(pC));
% 
% pC
