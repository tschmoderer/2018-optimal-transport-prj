clc
clear all
close all 

globals;

%% New implementation without staggered grid %%
N = 10; 
P = 10;
Q = 20; 

X = [0:N]/N; Y = [0:P]/P;
[XX,YY] = meshgrid(X,Y); YY = flipud(YY);

normalise = @(f) f/sum(f(:)); epsilon = 1e-10;
f0 = normalise(epsilon + gauss(0.3,0.3,0.05,0.05,N,P));
f1 = normalise(epsilon + gauss(0.8,0.8,0.05,0.05,N,P));

J = @(w) sum(sum(sum((w(:,:,:,1).^2 + w(:,:,:,2).^2)./max(w(:,:,:,3),max(epsilon,1e-10))))); % cost 

alpha = 1.0; % must be in ]0,2[
beta  = 1; % must be ine [0,1]
gamma = 2.0; % must be > 0

z  = zeros(P+1,N+1,Q,3);
w0 = zeros(P+1,N+1,Q,3); w1 = zeros(P+1,N+1,Q,3);

niter = 200;
cout = zeros(1,niter);
minF = zeros(1,niter);
divV = zeros(1,niter);

tic
for l = 1:niter
    w1 = w0 + alpha*(proxJ(2*z-w0,gamma) - z);
    [z, divV(l)] = projC(w1);
    w0 = w1;
    if mod(l,10) == 0
        contour(XX,YY,z(:,:,floor(Q/2)))
        title(['Iteration ',num2str(l)]);
        drawnow;
    end
    
    cout(l) = J(z);
    minF(l) = min(min(z(:,:,2)));
end
toc
close all

figure; 
surf(XX,YY,z(:,:,2));
title('Optimal transport');


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

