clear all
close all
clc
%% test new finite difference method
N = 10;
Q = 9;

X       = [0:N]/N; T = [0:Q]/Q;
[XX,YY] = meshgrid(X,T); YY = flipud(YY);

normalise = @(f) f/sum(f(:)); epsilon = 1e-10;
f0 = normalise(epsilon + gauss(0.5,0.05,N));
f1 = normalise(epsilon + gauss(0.5,0.05,N));

w = rand(Q+1,N+1,2);

dx  = @(m)  0.5*N*[zeros(Q+1,1), m(:,3:end) - m(:,1:end-2), zeros(Q+1,1)]; % dérivation centrée selon ---> x
dxS = @(dm) 0.5*N*[-dm(:,2), -dm(:,3), dm(:,2:end-3) - dm(:,4:end-1), dm(:,end-2), dm(:,end-1)];

% dx  = @(m) N*(m(:,[2:end end]) - m); % dérivation selon ---> x
dt  = @(f) Q*(f([2:end end],:) - f); % dérivation selon ---> t 
% dxS = @(dm) N*[-dm(:,1) , dm(:,1:end-2) - dm(:,2:end-1) , dm(:,end-1)];
dtS = @(df) Q*[-df(1,:) ; df(1:end-2,:) - df(2:end-1,:) ; df(end-1,:)];

%% check adjoint
rm = rand(Q+1,N+1); rdxrm = rand(Q+1,N+1);
dxrm = dx(rm); dxSrdxrm = dxS(rdxrm);

divM = sum(sum(sum(rm.*dxSrdxrm - dxrm.*rdxrm)))

rf = rand(Q+1,N+1); rdtrf = rand(Q+1,N+1);
dtrf = dt(rf); dtSrdtrf = dtS(rdtrf);

dtF = sum(sum(sum(rf.*dtSrdtrf - rdtrf.*dtrf)))
    
return
    %% construction opérateurs 
%     grad = @(w) [dx(w)]; 
     div  = @(w) -dxS(w(:,:,1));
    
     A    = @(w)  [div(w(:,:,1)) + dt(w(:,:,2)); w(end,:,2) ; w(1,:,2)];
     U    = @(y0,y1) [y1;zeros(Q-1,N+1);y0];
     AS   = @(Aw) cat(3,-dx(Aw(1:Q+1,:)),dtS(Aw(1:Q+1,:)) + U(Aw(Q+2,:),Aw(Q+3,:)));

       
    %% check adjoint
%     rw = rand(Q+1,N+1,2); rArw = rand(Q+3,N+1);
%     Arw = A(rw); ASrArw = AS(rArw);
% 
%     AAS = sum(sum(sum(rw.*ASrArw))) - sum(sum(sum(rArw.*Arw)))

    %% opérateurs de projetction
    y = [zeros(Q+1,N+1);f0;f1]; % second membre

    flat = @(x) x(:);
    resh = @(x) reshape(x,Q+3,N+1);
 
    do_cg =@(B,y) resh(cg(@(r)flat(B(resh(r))),y(:))); % solve B*x = y with CG
    pA = @(r) do_cg(@(s)A(AS(s)),r); % solve (A*A')*x = r

    pC = w + AS(pA(y - A(w)));
    
    err = @(w) norm(A(w)-y)/norm(y);
    error = err(pC);
    %% check div=0
%     w = rand(Q+1,N+1,2);
%     err = @(w) norm(A(w)-y)/norm(y);
%     fprintf('Error before projection: %.2e\n', err(w));
%     fprintf('Error before projection: %.2e\n', err(w + AS(pA(y - A(w)))));