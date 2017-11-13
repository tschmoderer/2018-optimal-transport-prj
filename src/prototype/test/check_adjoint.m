clc
clear all
close all

N = 10;
Q = 13;

% check interpolation

m = rand(Q+1,N+1);
f = rand(Q+1,N+1);
mbar = rand(Q+1,N+2);
fbar = rand(Q+2,N+1);

V = [reshape(m,(N+1)*(Q+1),1);reshape(f,(N+1)*(Q+1),1)];
U = [reshape(mbar,(N+2)*(Q+1),1);reshape(fbar,(N+1)*(Q+2),1)];

Vinterp = interpolation(U,N,Q);
Uinterpadj = interpolation_adj(V,N,Q);

e1 = (sum(U.*Uinterpadj) - sum(Vinterp.*V))/sum(Vinterp.*V)

% check divergence

d = rand(Q+1,N+1);
mbar = rand(Q+1,N+2);
fbar = rand(Q+2,N+1);

U = [reshape(mbar,(N+2)*(Q+1),1);reshape(fbar,(N+1)*(Q+2),1)];
d = reshape(d,(N+1)*(Q+1),1);

divU = divergence(U,N,Q);
divadjU = divergence_adj(d,N,Q);

e2 = (sum(divU.*d) - sum(divadjU.*U))/sum(divadjU.*U)

% check boundary

b = rand(2*(N+Q+2),1);
mbar = rand(Q+1,N+2);
fbar = rand(Q+2,N+1);

U = [reshape(mbar,(N+2)*(Q+1),1);reshape(fbar,(N+1)*(Q+2),1)];
boundU = boundary(U,N,Q);
boundUadj = boundary_adj(b,N,Q);

e3 = (sum(U.*boundUadj) - sum(U.*boundUadj))/sum(U.*boundUadj)







