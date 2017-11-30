clc
clear all
close all

N = 1;
Q = 1;

Gc = [0:N]/N;

m = zeros(N+1,Q+1);
f = zeros(N+1,Q+1);

mbar = zeros(N+2,Q+1);
fbar = zeros(N+1,Q+2);

%% interpolation
Im = (mbar(1:end-1,:) + mbar(2:end,:))/2;
If = (fbar(:,1:end-1) + fbar(:,2:end))/2;

%% divergence
dm = mbar(2:end,:) - mbar(1:end-1,:);
df = fbar(:,2:end) - fbar(:,1:end-1);
d = N*dm + Q*df;

%% boundary
bm1 = mbar(1,:); bm2 = mbar(end,:);
bf1 = fbar(:,1); bf2 = fbar(:,end);

b0m1 = zeros(1,Q+1); b0m2 = zeros(1,Q+1);
b0f1 = gauss(0.1,0.005,N)'; b0f2 = gauss(0.9,0.005,N)';

%% discrete functional
J = m.^2./f;
idx = find(m == 0 & f == 0);
J(idx) = 0;
J = sum(sum(J));

%% ProjC 


%% clear unusuful variables
clear idx