clc
clear all
close all 

globals;

%% New implementation without staggered grid %%
N = 48054; 
Q = 29; 

X       = [0:N]/N; T = [0:Q]/Q;
[XX,YY] = meshgrid(X,T); YY = flipud(YY);

normalise = @(f) f/sum(f(:)); epsilon = 1e-10;
f0 = normalise(epsilon + gauss(0.3,0.05,N));
f1 = normalise(epsilon + gauss(0.8,0.05,N));% + gauss(0.7,0.05,N));
%f1 = normalise(epsilon + 1./(1+10000*(X-0.5).^2));

%f0 = normalise(epsilon + indicatrix(0.2,0.8,N));
%f1 = normalise(epsilon + indicatrix(0.7,0.8,N) + indicatrix(0.3,0.5,N));

%!!% test musicale
f0 = audioread('f0.wav'); miniF0 = max(abs(f0));
f0 = normalise(epsilon + f0 + miniF0)';
f1 = audioread('f1.wav'); miniF1 = max(abs(f1));
f1 = normalise(epsilon + f1 + miniF1)';
%!!%

J = @(w) sum(sum(sum(w(:,:,1).^2./max(w(:,:,2),max(epsilon,1e-10))))); % cost 

alpha = 1.0; g = 2.0;

z  = zeros(Q+1,N+1,2);
w0 = zeros(Q+1,N+1,2); w1 = zeros(Q+1,N+1,2);
% t  = [Q:-1:0]/Q;
% tt = repmat(t',1,N+1);
% w0 = (1-tt).*repmat(f0,Q+1,1) + tt.*repmat(f1,Q+1,1);

niter = 2;
cout = zeros(1,niter);
minF = zeros(1,niter);
divV = zeros(1,niter);

tic
for l = 1:niter
    w1 = w0 + alpha*(proxJ(2*z-w0,g) - z);
    [z, divV(l)] = projC(w1);
    w0 = w1;
    if mod(l,10) == 0
      % surf(XX,YY,z(:,:,2))
        contour(XX,YY,z(:,:,2),35)
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

