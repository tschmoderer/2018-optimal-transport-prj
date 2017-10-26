function [m f] = projC(mbar,fbar)
global N;
global Q;
y = dlmread('files/Y.txt'); % size NxQ
y = reshape(y,Q+1,N+1);

[mybar fybar] = divergence_adjoint(y);
[my fy] = boundary_adjoint(y(:,1),y(:,end),y(1,:),y(end,:),N,Q);
end