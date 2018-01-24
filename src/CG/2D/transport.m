clc
clear all
close all 

globals;

%% New implementation without staggered grid %%
N = 21; 
P = 19;
Q = 20; 

X = (0:N)/N; Y = (0:P)/P; T = (0:Q)/Q;
[XX,YY] = meshgrid(X,Y);

normalise = @(f) f/sum(f(:)); epsilon = 0.05;
f0 = normalise(epsilon + gauss(0.2,0.2,0.1,N,P));
f1 = normalise(epsilon + gauss(0.8,0.8,0.1,N,P) + 0.6*gauss(0.7,0.4,0.07,N,P) + gauss(0.1,0.9,0.005,N,P));


% f0 = normalise(epsilon + indicatrix(0.2,0.8,0.1,0.16,N,P));
% f1 = normalise(epsilon + indicatrix(0.8,0.9,0.8,0.9,N,P));

alpha = 1.0; % must be in ]0,2[
beta  = 1.0; % must be ine [0,1] : 0 = interpolation, 1 = full transport
gamma = 0.5; % must be > 0

J = @(w) sum(sum(sum((w(:,:,:,1).^2 + w(:,:,:,2).^2)./(2*max(w(:,:,:,3),max(epsilon,1e-10))).^(beta)))); % cost 

z  = zeros(P+1,N+1,Q,3);
w0 = zeros(P+1,N+1,Q,3); w1 = zeros(P+1,N+1,Q,3);

niter = 1000;
cout = zeros(1,niter);
minF = zeros(1,niter);
divV = zeros(1,niter);

tic
for l = 1:niter
    w1 = w0 + alpha*(proxJ(2*z-w0,gamma,beta) - z);
    [z, divV(l)] = projC(w1);
    w0 = w1;
    
    cout(l) = J(z);
    minF(l) = min(min(z(:,:,2)));    
    
    % Affichage
    if mod(l,50) == 0
        contour(XX,YY,z(:,:,floor(4*Q/5)))
        title(['Iteration ',num2str(l)]);
        drawnow;
    end
end
toc

close all
I = floor(Q/6);
figure; 
subplot(2,3,1), 
contour(XX,YY,z(:,:,1,3));
title('initial density - t=0');

for i = 2:5
   subplot(2,3,i),
   contour(XX,YY,z(:,:,i*I,3));
   title(['density at t = ',num2str(T(i*I))]);
end

subplot(2,3,6), 
contour(XX,YY,z(:,:,end,3));
title('Target density - t=1');

figure; 
subplot(311)
plot([1:niter],cout);
title('cout');
subplot(312)
plot([1:niter],minF);
title('Minimum de F');
subplot(313);
plot([1:niter],divV);
title('divergence violation');

