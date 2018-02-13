clc
clear all
close all 

globals;

%% New implementation without staggered grid %%
N = 20; 
P = 20;
Q = 20; 

X = (0:N)/N; Y = (0:P)/P; T = (0:Q)/Q;
[XX,YY] = meshgrid(X,Y);

normalise = @(f) f/sum(f(:)); epsilon = 0.05;
f0 = normalise(epsilon + gauss(0.2,0.2,0.05,N,P));
f1 = normalise(epsilon + gauss(0.8,0.8,0.05,N,P));
eps = min(f0(:));

alpha = 1.0; % must be in ]0,2[
beta  = 1.0; % must be ine [0,1] : 0 = interpolation, 1 = full transport
gamma = 1.0; % must be > 0

J = @(w) sum(sum(sum((w(:,:,:,1).^2 + w(:,:,:,2).^2)./(2*max(w(:,:,:,3),max(epsilon,1e-10))).^(beta)))); % cost 

zV  = zeros(P+1,N+1,Q,3); wV0 = zeros(P+1,N+1,Q,3); wV1 = zeros(P+1,N+1,Q,3);
zU  = zeros(P+2,N+2,Q+1,3); wU0 = zeros(P+2,N+2,Q+1,3); wU1 = zeros(P+2,N+2,Q+1,3);

T = ([Q-1:-1:-1]+0.5)/(Q-1); T = ([-1:Q-1]+0.5)/(Q-1); TT = zeros(P+2,N+2,Q+1); F0 = []; F1 = []; 
for i = 1:Q+1, TT(:,:,i) = T(i)*ones(P+2,N+2); F0 = cat(3,F0,[f0 ones(P+1,1); ones(1,N+2)]); F1 = cat(3,F1,[f1 ones(P+1,1); ones(1,N+2)]);end
wU0(:,:,:,3) = (1-TT).*F0 + TT.*F1;
wV0 = interp(wU0); 
zU = wU0; zV = wV0;

niter = 1000;
cout = zeros(1,niter);
minF = zeros(1,niter);
divV = zeros(1,niter);

tic
for l = 1:niter	
    wU1 = wU0 + alpha*(projC(2*zU - wU0) - zU);
    wV1 = wV0 + alpha*(proxJ(2*zV - wV0,gamma,beta) - zV);
    zU  = projCs(wU1,wV1);
    zV  = interp(zU);
    
    wU0 = wU1;
    wV0 = wV1;
    
    cout(l) = J(zV);
    minF(l) = min(min(min(zV(:,:,:,3))));    
    
    % Affichage
    if mod(l,20) == 0
        contour(XX,YY,zV(:,:,floor(4*Q/5),3))
        title(['Iteration ',num2str(l)]);
        drawnow;
    end
end
toc

close all
I = floor(Q/6);
figure; 
subplot(2,3,1), 
contour(XX,YY,zV(:,:,1,3));
title('initial density - t=0');

for i = 2:5
   subplot(2,3,i),
   contour(XX,YY,zV(:,:,i*I,3));
   title(['density at t = ',num2str(T(i*I))]);
end

subplot(2,3,6), 
contour(XX,YY,zV(:,:,end,3));
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







% Affichage images
I = floor(Q/16);
figure; 
subplot(4,4,1), 
imshow(mat2gray(zV(:,:,1,3)));
title('initial density - t=0');

for i = 2:15
   subplot(4,4,i),
   imshow(mat2gray(zV(:,:,i*I,3)));
   title(['density at t = ',num2str(T(i*I))]);
end

subplot(4,4,16), 
imshow(mat2gray(zV(:,:,end,3)));
title('Target density - t=1');




