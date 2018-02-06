clc
clear all
close all 

% Programme principal
% On essaye d'agir que sur l'intérieur de la grille 
% Timothée Schmoderer 
% INSA Rouen Normandie 2017/2018

globals;

%% New implementation without staggered grid %%
N = 21; 
Q = 19; 

X       = (0:N)/N; T = (0:Q)/Q;
[XX,YY] = meshgrid(X,T); YY = flipud(YY);

normalise = @(f) f/sum(f(:)); epsilon = 1e-10;
obstacle = zeros(Q+1,N+1);

f0 = normalise(epsilon + gauss(0.2,0.05,N));
f1 = normalise(epsilon + gauss(0.8,0.05,N));
epsilon = min(f0);

alpha = 1.0; % must be in ]0,2[
beta  = 1.0; % must be ine [0,1]
gamma = 1.0; % must be > 0

J = @(w) sum(sum(sum(w(:,:,1).^2./max(w(:,:,2),max(epsilon,1e-10)).^beta))); % cost 

z  = zeros(Q+1,N+1,2);
t  = (Q:-1:0)/Q; tt = repmat(t',1,N+1);
w0 = (1-tt).*repmat(f0,Q+1,1) + tt.*repmat(f1,Q+1,1);
w1 = w0;

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


