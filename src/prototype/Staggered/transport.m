clc
clear all
close all 

globals;

N = 20; Q = 21; niter = 1000;
eps = 1e-10; alpha = 1.0; g = 1.0; b = 1.0;

f0 = zeros(1,N+1); f1 = zeros(1,N+1);
zV = zeros(Q+1,N+1,2); wV0 = zeros(Q+1,N+1,2); wV1 = zeros(Q+1,N+1,2);
zU = zeros(Q+2,N+2,2); wU0 = zeros(Q+2,N+2,2); wU1 = zeros(Q+2,N+2,2);
obstacle = zeros(Q+1,N+1);
cout = zeros(1,niter); minF = zeros(1,niter);

J = @(w) 0.5*sum(sum(sum(w(:,:,1).^2./max(w(:,:,2),max(eps,1e-10)).^b)));

normalise = @(f) f/sum(f(:));
	 
f0 = normalise(eps + gauss(0.2,0.05));
f1 = normalise(eps + gauss(0.8,0.05));

eps = min(f0(:));

for i = 1:niter
    wU1 = wU0 + alpha*(projC(2*zU - wU0) - zU);
    wV1 = wV0 + alpha*(proxJ(2*zV - wV0,b,g,obstacle) - zV);
    zU  = projCs(wU1,wV1);
    zV  = interp(zU);

    wU0 = wU1;
    wV0 = wV1;

    cout(i) = J(zV);
    minF(i) = min(min(zV(:,:,2)));
    if (mod(i,10) == 0) 
        fprintf('Iteration : %d, cout : %f\n',i, cout(i));
    end  
end
