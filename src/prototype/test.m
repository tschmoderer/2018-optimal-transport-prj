% check value given by the fortran code
close all

N=30;Q=30;
%m = reshape(dlmread('Fortran/solution/m'),Q+1,N+1);
%f = reshape(dlmread('Fortran/solution/f'),N+1,Q+1)';
%m = reshape(dlmread('m'),Q+1,N+1)';
%f = reshape(dlmread('f'),N+1,Q+1)';
%mbar = reshape(dlmread('mbar'),Q+1,N+2);
%fbar = reshape(dlmread('fbar'),N+1,Q+2)';

[X,Y] = meshgrid(linspace(0,1,N+1),linspace(0,1,Q+1));
%surf(X,Y,f);
xlabel('x')
ylabel('t')
zlabel('f/m')
%test svg Y
YY = reshape(dlmread('Fortran/files/Y'),N+1,Q+1)';
ZZ = reshape(dlmread('Fortran/files/solution'),N+1,Q+1)';
%proxJ(m,f,gamma);
surf(X,Y,YY)
%surf(X,Y,ZZ)
xlabel('x')
ylabel('t')