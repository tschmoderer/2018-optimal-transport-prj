clc
clear all
close all

N = 2; P = 3; Q = 2;

%% Centered Data
% Space time cube : 

GcX = (0:N)/N;  GcY = (0:P)/P; GcT = (0:Q)/Q;
Gc = meshgrid(GcX,GcY,GcT);
% y
% ^
% |
% | (P+1) pts
% |
%  --------- > x
%   (N+1) pts

% Variable
m = rand([size(Gc) 2]); 
f = rand(size(Gc));
V = cat(4,m,f); % V(:,:,:,1:2) = m, V(:,:,:,3) = f

%% Staggered Data
% 3 space time cube
GsX = meshgrid(((-1:N)+0.5)/N,GcY,GcT); 
GsY = meshgrid(GcX,((-1:P)+0.5)/P,GcT);
GsT = meshgrid(GcX,GcY,((-1:Q)+0.5)/Q);

% Variable
mbar1 = rand(size(GsX)); mbar2 = rand(size(GsY));
fbar = rand(size(GsT));

U.mbar1 = mbar1; U.mbar2 = mbar2; U.fbar = fbar;

%% Opérateurs
% Interpolation 

I = @(U) 0.5*cat(4,0.5*(U.mbar1(1:end-1,:,:) - U.mbar1(2:end,:,:)) + 0.5*(U.mbar2(:,1:end-1,:) - U.mbar2(:,2:end,:)),U.fbar(:,:,1:end-1)-U.fbar(:,:,2:end));

% divergence

div = @(mbar,fbar) N*(mbar(2:end,:,:,1) - mbar(1:end-1,:,:,1)) + P*(mbar(:,2:end,:,2) - mbar(:,1:end-1,:,2)) + Q*(fbar(:,:,2:end) - fbar(:,:,1:end-1));


