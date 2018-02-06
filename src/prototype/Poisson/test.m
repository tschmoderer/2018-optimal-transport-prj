clear all
close all
clc

%% test DF
N = 19;
Q = 21;

normalise = @(f) f/sum(f(:)); epsilon = 1e-10;
f0 = normalise(epsilon + gauss(0.5,0.05,N));
f1 = normalise(epsilon + gauss(0.5,0.05,N));

dx  = @(m) N*(m(:,[2:end end]) - m); % dérivation selon ---> x
dt  = @(f) Q*(f([2:end end],:) - f); % dérivation selon ---> t 
dxS = @(dm) N*[-dm(:,1) , dm(:,1:end-2) - dm(:,2:end-1) , dm(:,end-1)];
dtS = @(df) Q*[-df(1,:) ; df(1:end-2,:) - df(2:end-1,:) ; df(end-1,:)];

div  = @(w) -dxS(w(:,:,1)) + dt(w(:,:,2));
grad = @(d) cat(3,-dx(d(1:Q+1,:)),dtS(d(1:Q+1,:)));

% A  = @(w)  [-dxS(w(:,:,1)) + dt(w(:,:,2)); w(end,:,2); w(1,:,2)];
% AS = @(Aw) cat(3,-dx(Aw(1:Q+1,:)),dtS(Aw(1:Q+1,:)) + [Aw(Q+3,:);zeros(Q-1,N+1);Aw(Q+2,:)]);

y = [zeros(Q+1,N+1);f0;f1];

w = zeros(Q+1,N+1,2); tmpw = w;

dw = div(w);
rdw = rand(size(dw)); dSrdw = grad(rdw);

testDiv = sum(sum(sum(w.*dSrdw))) - sum(sum(dw.*rdw))


% find u such that -DELTA(u) = -div(w) := b
w(:,[1 end],1) = 0; w([1 end],:,2) = [f1 ; f0];
b = div(w);

u = poisson_neumann_2d(-b,f0,f1,N,Q);

% gu = grad(u); gu0 = gu; gu0(2:end-1,2:end-1) = 0;
% pCw = w + gu - gu0


u2 = poisson2d_Neumann(-b); 
gu2 = grad(u2); gu02 = gu2; gu02(2:end-1,2:end-1) = 0;
pCw2 = w + gu2 - gu02

figure; 
subplot(221) 
surf(pCw(:,:,1)); title('pCw - m'); 
subplot(222) 
surf(pCw(:,:,2)); title('pCw - f'); 

subplot(223) 
surf(pCw2(:,:,1)); title('pCw2 - m'); 
subplot(224) 
surf(pCw2(:,:,2)); title('pCw2 - f'); 