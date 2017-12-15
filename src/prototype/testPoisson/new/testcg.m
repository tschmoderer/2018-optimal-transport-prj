clc
clear all
close all 

%% test Conjugate gradient
A = @(x) 2*x;
mustBe50 = cg(A,100)
A = @(x) [x(1)+2*x(2);2*x(2)+x(1)];
mustBe1thirdtwice = cg(A,[1;1])
mustBe0twice = cg(A,[0;0])

% couple = cg(A,[1 0;1 0])

% return
%% use conjugate gradient
globals;

N = 12;
Q = 13;

w = rand(Q+1,N+1,2);
normalise = @(f) f/sum(f(:)); epsilon = 1e-10;
f0 = normalise(epsilon + gauss(0.3,0.05,N));
f1 = normalise(epsilon + gauss(0.8,0.05,N));

dx  = @(m) (m(:,[2:end end]) - m); % dérivation selon ---> x
dt  = @(f) (f([2:end end],:) - f); % dérivation selon ---> t 
dxS = @(dm) [-dm(:,1) , dm(:,1:end-2) - dm(:,2:end-1) , dm(:,end-1)];
dtS = @(df) [-df(1,:) ; df(1:end-2,:) - df(2:end-1,:) ; df(end-1,:)];

%% construction opérateurs 
 div  = @(w) -dxS(w(:,:,1));

 A    = @(w)  [div(w(:,:,1)) + dt(w(:,:,2)); w(end,:,2) ; w(1,:,2)];
 U    = @(y0,y1) [y1;zeros(Q-1,N+1);y0];
 AS   = @(Aw) cat(3,-dx(Aw(1:Q+1,:)),dtS(Aw(1:Q+1,:)) + U(Aw(Q+2,:),Aw(Q+3,:)));

%% opérateurs de projetction

y = [zeros(Q+1,N+1);f0;f1]; % second membre
opts.epsilon = 1e-9;
opts.niter_max = 150;
flat = @(x) x(:);
resh = @(x) reshape(x,Q+3,N+1);
%cg = @(B,y) resh(perform_cg(@(r)flat(B(resh(r))),y(:),opts)); % ma fonction gradient conjugué

pA = @(r) cg(@(s)A(AS(s)),r);
pA = @(r) resh(cg(@(r)flat(A(AS(resh(r)))),y(:)))
pC = w + AS(pA(y - A(w)));

err = @(w) norm(A(w)-y)/norm(y);
error = err(pC)




%% check div=0
%     w = rand(Q+1,N+1,2);
%     err = @(w) norm(A(w)-y)/norm(y);
%     fprintf('Error before projection: %.2e\n', err(w));
%     fprintf('Error before projection: %.2e\n', err(w + AS(pA(y - A(w)))));