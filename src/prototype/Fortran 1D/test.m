clear all 
close all
clc

%% Gaussienne 
gauss = @(mu,sigma,N) exp(-0.5*(([0:N]/N - mu)/sigma).^2);

%% Normalise 
normalise = @(f) f/sum(f(:));

%% Dérivation selon x
dx = @(m,N,Q) N*[m(:,2:N+1) - m(:,1:N), zeros(Q+1,1)];	

%%Dérivation selon t
dt = @(f,N,Q) Q*[f(2:Q+1,:) - f(1:Q,:) ; zeros(1,N+1)];

%% Adjoint de dx
dxS = @(dm,N,Q) N*[-dm(:,1), dm(:,1:N-1) - dm(:,2:N), dm(:,N)];

%% Adjoint de dt 
dtS = @(df,N,Q) Q*[-df(1,:); df(1:Q-1,:) - df(2:Q,:); df(Q,:)];

%% Opérateur A
A = @(w,N,Q) [-dxS(w(:,:,1),N,Q) + dt(w(:,:,2),N,Q); w(Q+1,:,2); w(1,:,2)];

%% Opérateur adjoint de A
AS = @(Aw,N,Q) cat(3,-dx(Aw(1:Q+1,:),N,Q),dtS(Aw(1:Q+1,:),N,Q) + [Aw(Q+3,:);zeros(Q-1,N+1);Aw(Q+2,:)]);
	
%% Function aplatir
flat = @(x,N,Q) reshape(x,(Q+3)*(N+1),1);

%% Fonction reshape
resh = @(x,N,Q) reshape(x,Q+3,N+1);	

%% Projection sur C
projC = @(w,y,N,Q) w + AS(resh(cg(@(s,N,Q) flat(A(AS(resh(s,N,Q),N,Q),N,Q),N,Q), flat(y-A(w,N,Q),N,Q),N,Q),N,Q),N,Q); 

%% le cout
% 	function cost(w,N,Q) result(c)
% 		integer :: N,Q
% 		double precision, dimension(1:Q+1,1:N+1,2) :: w
% 		double precision:: c
% 		c = sum(w(:,:,1)**2/w(:,:,2))
% 	end function cost


%% TEST %%
% N = 30; Q = 30;
% rm = rand(Q+1,N+1); rdxrm = rand(Q+1,N+1);
% rf = rand(Q+1,N+1); rdtrf = rand(Q+1,N+1);
% 
% dxrm     = dx(rm,N,Q);
% dxSrdxrm = dxS(rdxrm,N,Q);
% dtrf     = dt(rf,N,Q);
% dtSrdtrf = dtS(rdtrf,N,Q);
% 
% divM = sum(sum(rm.*dxSrdxrm - dxrm.*rdxrm))
% dT   = sum(sum(rf.*dtSrdtrf - dtrf.*rdtrf))
% 
% rw = rand(Q+1,N+1,2); rArw = rand(Q+3,N+1);
% 
% Arw    = A(rw,N,Q);
% ASrArw = AS(rArw,N,Q);
% 
% AAS = sum(sum(sum(rw.*ASrArw))) - sum(sum(rArw.*Arw)) 
% 
% 	
% y = zeros(Q+3,N+1);
% % 	!! test projection
% epsilon = 1e-10;
% f0 = normalise(epsilon + gauss(0.5d0,0.05d0,N));
% f1 = normalise(epsilon + gauss(0.5d0,0.05d0,N));
% y(1:Q+1,:) = 0;
% y(Q+2,:) = f0;
% y(Q+3,:) = f1;
% 
% w = rand(Q+1,N+1,2);
% BP = sum(sum((A(w,N,Q) - y).^2))/sum(sum(y.^2))
% pC = w + AS(resh(cg(@(s,N,Q) flat(A(AS(resh(s,N,Q),N,Q),N,Q),N,Q), flat(y-A(w,N,Q),N,Q),N,Q),N,Q),N,Q);
% AP = sum(sum((A(pC,N,Q) - y).^2))/sum(sum(y.^2))

N = 21;
Q = 19;
epsilon = 1e-10; alpha = 1.0; gamma = 1.0;
w0 = zeros(Q+1,N+1,2); w1 = zeros(Q+1,N+1,2); z = zeros(Q+1,N+1,2);

f0 = normalise(epsilon + gauss(0.5d0,0.05d0,N));
f1 = normalise(epsilon + gauss(0.5d0,0.05d0,N));

y = zeros(Q+3,N+1);
y(1:Q+1,:) = 0;
y(Q+2,:) = f0;
y(Q+3,:) = f1;

for i = 1:2000
    w1 = w0 + alpha*(proxJ(2*z-w0,gamma,N,Q) - z);
    z = projC(w1,y,N,Q);
    if mod(i,10) == 0
        surf(z(:,:,2),'Edgecolor','none')
        title(['Iteration : ',num2str(i)])
        drawnow
    end
    w0 = w1;
end