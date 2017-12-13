clc
clear all
close all 

%% New implementation without staggered grid %%
N = 10; 
Q = 9; 

m = zeros(Q+1,N+1);
f = zeros(Q+1,N+1); 

w = zeros(Q+1,N+1,2);

w(:,:,1:end-1) = m;
w(:,:,end)     = f;

normalise = @(f) f/sum(f(:));
f0 = normalise(0.05 + gauss(0.2,0.005,N));
f1 = normalise(0.05 + gauss(0.8,0.005,N));

J = @(w) w(:,:,1).^2./w(:,:,2); % /!\ works only in 2d


alpha = 1.998; g = 1.0;

z0 = zeros(Q+1,N+1,2); z1 = zeros(Q+1,N+1,2);
w0 = zeros(Q+1,N+1,2); w1 = zeros(Q+1,N+1,2);

for l = 1:10
    w1 = w0 + alpha*(proxJ(2*z0-w0,g) - z0);
    z1 = projC(w1,f0,f1,N,Q);
    w0 = w1;
    z0 = z1;
    
    surf(z0(:,:,2))
    title(['Iteration ',num2str(l)]);
    drawnow;
end
