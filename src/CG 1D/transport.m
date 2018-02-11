clc
clear all
close all 

% Programme principal
% Algorithme de Douglas - Rachford pour le trans port optimal 
% Timothée Schmoderer 
% INSA Rouen Normandie 2017/2018

globals;

%% New implementation without staggered grid %%
N = 100; 
Q = 100; 

X       = (0:N)/N; T = (0:Q)/Q;
[XX,YY] = meshgrid(X,T); YY = flipud(YY);

normalise = @(f) f/sum(f(:)); epsilon = 1e-10;
obstacle = zeros(Q+1,N+1);

%% Test 1 %%
f0 = normalise(epsilon + gauss(0.2,0.05,N));
f1 = normalise(epsilon + gauss(0.8,0.05,N));

%% Test 1.5 %%
% X0 = [0.0 0.1 0.2 0.4 0.7 1.0]; Y0 = [0.1 1.5 0.3 0.4 0.5 0.1];
% X = linspace(0,1,N+1);
% f0 = spline(X0,Y0,X); f0 = normalise(eps + f0 - min(f0));
% 
% X1 = [0.0 0.1 0.2 0.4 0.7 1.0]; Y1 = [0.1 0.5 0.3 0.4 1.5 0.1];
% f1 = spline(X1,Y1,X); f1 = normalise(eps + f1 - min(f1));

% %% Test 2 %%
% f0 = normalise(epsilon + gauss(0.2,0.05,N));
% f1 = normalise(epsilon + gauss(0.8,0.05,N) + gauss(0.5,0.05,N));
% 
% %% Test 3 %%
% f0 = normalise(epsilon + indicatrix(0.2,0.8,N));
% f1 = normalise(epsilon + indicatrix(0.7,0.8,N) + indicatrix(0.3,0.5,N));
% 
% %% Test 4 %%
% f0 = normalise(epsilon + gauss(0.2,0.05,N));
% f1 = normalise(epsilon + gauss(0.8,0.05,N));

obstacle = zeros(Q+1,N+1); % 0 : no obstacle, 1 : obstacle
% obstacle(15,1:end) = 1; obstacle(15,4:6) = 0;
%obstacle(20,1:end) = 1; obstacle(20,16:18) = 0;
% obstacle(7,1:end)  = 1; obstacle(7,12:14) = 0;
% 
% %% Test 5 %%
% obstacle(45,1:end) = 1; obstacle(45,40:45) = 0;
% obstacle(30,1:end) = 1; obstacle(30,4:6) = 0;

J = @(w) sum(sum(sum(w(:,:,1).^2./w(:,:,2)))); % cost 

alpha = 1.0; % must be in ]0,2[
beta  = 1.0; % must be ine [0,1]
gamma = 1.0; % must be > 0

z  = zeros(Q+1,N+1,2);
w0 = zeros(Q+1,N+1,2); w1 = zeros(Q+1,N+1,2);
% t  = (Q:-1:0)/Q;
% tt = repmat(t',1,N+1);
% w0 = (1-tt).*repmat(f0,Q+1,1) + tt.*repmat(f1,Q+1,1);

niter = 200;
cout = zeros(1,niter);
minF = zeros(1,niter);
divV = zeros(1,niter);

tic
for l = 1:niter
    w1 = w0 + alpha*(proxJ(2*z-w0,beta,gamma,obstacle) - z);
    [z, divV(l)] = projC(w1);

    w0 = w1;
    
    cout(l) = J(z);
    minF(l) = min(min(z(:,:,2)));

    % Affichage
    if mod(l,20) == 0
        contour(XX,YY,z(:,:,2),35)
        title(['Iteration ',num2str(l)]);
        drawnow;
    end
end
toc
close all

figure; 
surf(XX,YY,z(:,:,2),'EdgeColor','none');
title('Transport optimal');
xlabel('x')
ylabel('t')


figure; 
plot(X,f0)
title('Densité initiale');
xlabel('x')

figure; 
plot(X,f1)
title('Densité cible');
xlabel('x')

figure; 
contour(XX,YY,z(:,:,2),50);
title('Transport optimal');
xlabel('x')
ylabel('t')
if sum(sum(obstacle)) > 0
 hold on; 
 contour(XX,YY,obstacle,50,'LineColor','k');
end

figure; 
subplot(311)
plot((1:niter),cout);
title('Energie');
xlabel('Iteration')
subplot(312)
plot([1:niter],minF);
title('Minimum de la densité');
subplot(313);
plot([1:niter],divV);
title('Violation de la contrainte div = 0');


