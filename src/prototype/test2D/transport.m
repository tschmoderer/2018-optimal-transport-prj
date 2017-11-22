clc
clear all
close all

globals; 

% %% Initialisation %%

N = 2;
P = 3;
Q = 2;

 divergence;
 boundary;
 interpolation;
% 
% matrices projection sur C %
 A = [D ; B]; 
 delta = A*A'; 
 
 normalise = @(f) f/sum(f(:));
 
 sigma = 0.05; mini = 0.0001;
 f0 = normalise(gauss(0.5,0.5,sigma,sigma,N,P,mini));
 f1 = normalise(gauss(0.5,0.5,sigma,sigma,N,P,mini));

y = [zeros((N+1)*(P+1)*(Q+1),1) ; zeros(2*2*(N+P+2)*(Q+1),1) ; reshape(f0,[],1) ; reshape(f1,[],1)];
Cst = A'*(delta\y);
 
PC = eye(2*(N+2)*(P+2)*(Q+1) + (N+1)*(P+1)*(Q+2)) - A'*(delta\A);

% matrice de projection sur G2
pG2 = inv(eye(2*(N+2)*(P+2)*(Q+1) + (N+1)*(P+1)*(Q+2)) + Interp'*Interp);

m = zeros(N+1,P+1,Q+1,2);
f = zeros(N+1,P+1,Q+1);

V = [reshape(m(:,:,:,1),[],1) ; reshape(m(:,:,:,2),[],1) ; reshape(f,[],1)];

mbar = zeros(N+2,P+2,Q+1,2);
fbar = zeros(N+1,P+1,Q+2);

U = [reshape(mbar(:,:,:,1),[],1) ; reshape(mbar(:,:,:,2),[],1) ; reshape(fbar,[],1)];


