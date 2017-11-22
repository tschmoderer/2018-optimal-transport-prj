clc
clear all
close

%% Test Divergence %%
N = 1;
P = 2;
Q = 1;

mbar = rand(P+2,N+2,Q+1,2);
fbar = rand(P+1,N+1,Q+2);

U = [reshape(mbar(:,:,:,1),[],1) ; reshape(mbar(:,:,:,2),[],1) ; reshape(fbar,[],1)];

divergence;

% reference %
dX = N*(mbar(:,2:end,:,1) - mbar(:,1:end-1,:,1)); % div selon x
dY = P*(mbar(2:end,:,:,2) - mbar(1:end-1,:,:,2)); % div selon y
d = dX(2:end,:,:) + dX(1:end-1,:,:) + dY(:,2:end,:) + dY(:,1:end-1,:) + Q*(fbar(:,:,2:end) - fbar(:,:,1:end-1)); % div totale


e1 = sum(sum(sum(reshape(D*U,[P+1,N+1,Q+1]) - d)))


%% test Interpolation %%
interpolation;

% reference %
InterpM1 = mbar(:,:,:,1) + mbar(:,:,:,1);
InterpF = fbar(:,:,2:end) + fbar(:,:,1:end-1);

% a revoir
%e2 = sum(sum(sum(reshape(I*U,[P+1,N+1,Q+1,2]) - Interp)))


%% test boundary %%

% mbar(:,:,1,1) = [3 6 8;3 5 6;7 7 6];
% mbar(:,:,2,1) = [6 8 10;6 10 4;1 4 9];
% mbar(:,:,1,2) = [3 6 8;3 5 6;7 7 6];
% mbar(:,:,2,2) = [6 8 10;6 10 4;1 4 9];
% 
% 
% fbar (:,:,1) = [1 2 ; 3 4];
% fbar (:,:,2) = [5 6 ; 7 8];
% fbar (:,:,3) = [9 10 ; 11 12];
% 
% U = [reshape(mbar(:,:,:,1),[],1) ; reshape(mbar(:,:,:,2),[],1) ; reshape(fbar,[],1)];

boundary;

B*U


%% test contraintes sur C %%
A = [D ; B];
delta = A*A';

% test cout 
m = rand(P+1,N+1,Q+1,2);
f = rand(P+1,N+1,Q+1);

V = [reshape(m(:,:,:,1),[],1) ; reshape(m(:,:,:,2),[],1) ; reshape(f,[],1)];

R = (V(1:(N+1)*(P+1)*(Q+1)).^2 + V((N+1)*(P+1)*(Q+1)+1:2*(N+1)*(P+1)*(Q+1)).^2);
R = R./V(2*(N+1)*(P+1)*(Q+1)+1:end);
R = sum(R);

%%
























