clc
clear all
%close all

globals;

N = 7; Q = 5; niter = 200;
eps = 1e-10; alpha = 1.0; g = 1.0; b = 1;

normalise = @(f) f/sum(f(:));

f0 = normalise(eps + gauss(0.2,0.05,N));
f1 = normalise(eps + gauss(0.8,0.05,N));

zV = zeros(Q+1,N+1,2); wV0 = zeros(Q+1,N+1,2); wV1 = zeros(Q+1,N+1,2);
zU = zeros(Q+2,N+2,2); wU0 = zeros(Q+2,N+2,2); wU1 = zeros(Q+2,N+2,2); 

J = @(w) 0.5*sum(sum(w(:,:,1).^2./max(w(:,:,2),max(eps,1e-10)).^b));

obstacle = zeros(Q+1,N+1);
cout = zeros(1,niter); minF = zeros(1,niter);

T = ([Q:-1:-1]+0.5)/Q; TT = repmat(T',1,N+1);
wU0(:,:,2) = [(1-TT).*repmat(f0,Q+2,1) + TT.*repmat(f1,Q+2,1) zeros(Q+2,1)];
wV0 = interp(wU0); 
zU = wU0; zV = wV0;

for i = 1:niter
    wU1 = wU0 + alpha*(projC(2*zU - wU0) - zU);
    wV1 = wV0 + alpha*(proxJ(2*zV - wV0,b,g,obstacle) - zV);
    zU  = projCs(wU1,wV1);
    zV  = interp(zU);

    wU0 = wU1;
    wV0 = wV1;

    cout(i) = J(zV);

    if (mod(i,1) == 0) 
        cout(i)
        surf(zV(:,:,2))
        title(['Iteration : ',num2str(i)])
        drawnow
    end
    minF(i) = min(min(zV(:,:,2)));
end