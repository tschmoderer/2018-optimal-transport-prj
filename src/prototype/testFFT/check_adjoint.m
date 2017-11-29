clc
clear all
close all

globals;

N = 2;
Q = 2;

% Matrice de l'opérateur b %
Bm = zeros(2*(Q+1),(N+2)*(Q+1));
Bm(1:Q+1,1:Q+1) = eye(Q+1);
Bm(Q+2:end,end-Q:end) = eye(Q+1);

Bf = [];
for i = 1:N+1
    Bf = blkdiag(Bf,[1 zeros(1,Q+1);zeros(1,Q+1) 1]);
end

B = blkdiag(Bm,Bf);
% matrice d'interpolation %

Im = zeros((N+1)*(Q+1),(N+2)*(Q+1));
for i = 1:(N+1)*(Q+1)
    for j = 1:(N+2)*(Q+1)
        if i == j 
            Im(i,j) = 1;
        elseif j == i+Q+1
            Im(i,j) = 1;
        end
    end
end
dia = zeros(Q+1,Q+2);
for i = 1:Q+1
    for j = 1:Q+2
        if i == j 
            dia(i,j) = 1;
        elseif j == i+1
            dia(i,j) = 1;
        end
    end
end
If = [];
for i = 1:N+1
    If = blkdiag(If,dia);
end

Interp = 0.5*blkdiag(Im,If);

% matrice de projection sur G2
pG2 = inv(eye((N+1)*(Q+2)+(N+2)*(Q+1)) + Interp'*Interp);

% Matrice de la divergence %
Dm = zeros((N+1)*(Q+1),(N+2)*(Q+1));
for i = 1:(N+1)*(Q+1)
    for j = 1:(N+2)*(Q+1)
        if i == j 
            Dm(i,j) = -1;
        elseif j == i+Q+1
            Dm(i,j) = 1;
        end
    end
end
dia = zeros(Q+1,Q+2);
for i = 1:Q+1
    for j = 1:Q+2
        if i == j 
            dia(i,j) = -1;
        elseif j == i+1
            dia(i,j) = 1;
        end
    end
end
Df = [];
for i = 1:N+1
    Df = blkdiag(Df,dia);
end

D = [N*Dm Q*Df];

% matrices projection sur C %
A = [D ; B]; 
delta = A*A'; 

sigma = 0.05; mini = 0.0001;
f0 = gauss(0.2,sigma,N,mini); 
f1 = gauss(0.8,sigma,N,mini); 

y = [zeros((N+1)*(Q+1),1) ; zeros(2*(Q+1),1) ; reshape([f1;f0],2*(N+1),1)];
Cst = A'*(delta\y);

P = eye((N+1)*(Q+2)+(N+2)*(Q+1)) - A'*(delta\A);

% check interpolation

m = rand(Q+1,N+1);
f = rand(Q+1,N+1);
mbar = rand(Q+1,N+2);
fbar = rand(Q+2,N+1);

V = [reshape(m,(N+1)*(Q+1),1);reshape(f,(N+1)*(Q+1),1)];
U = [reshape(mbar,(N+2)*(Q+1),1);reshape(fbar,(N+1)*(Q+2),1)];

% Vinterp = interpolation(U,N,Q);
% Uinterpadj = interpolation_adj(V,N,Q);
Vinterp = Interp*U;
Uinterpadj = Interp'*V;

e1 = (sum(U.*Uinterpadj) - sum(Vinterp.*V))/sum(Vinterp.*V)

% check divergence

d = rand(Q+1,N+1);
mbar = rand(Q+1,N+2);
fbar = rand(Q+2,N+1);

U = [reshape(mbar,(N+2)*(Q+1),1);reshape(fbar,(N+1)*(Q+2),1)];
d = reshape(d,(N+1)*(Q+1),1);

% divU = divergence(U);
% divadjU = divergence_adj(d);

divU = D*U;
divadjU = D'*d;

e2 = (sum(divU.*d) - sum(divadjU.*U))/sum(divadjU.*U)

% check boundary

b = rand(2*(N+Q+2),1);
mbar = rand(Q+1,N+2);
fbar = rand(Q+2,N+1);

U = [reshape(mbar,(N+2)*(Q+1),1);reshape(fbar,(N+1)*(Q+2),1)];
% boundU = boundary(U);
% boundUadj = boundary_adj(b);
boundU = B*U;
boundUadj = B'*b;

e3 = (sum(U.*boundUadj) - sum(U.*boundUadj))/sum(U.*boundUadj)


% check opérateurs 

N = 2; Q = 2;

mbar = [1 9 10 2;2 3 4 8; 4 10 10 6];
fbar = [7 7 10;5 9 6;3 7 3;2 5 7];

U = [reshape(mbar,(N+2)*(Q+1),1);reshape(fbar,(N+1)*(Q+2),1)];

D*U;
