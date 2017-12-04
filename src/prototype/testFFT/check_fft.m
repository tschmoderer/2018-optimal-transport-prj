clc
close all

N = 5; Q = 5; P = Q;

%% reference
% Matrice de l'opérateur b %
Bm = zeros(2*(Q+1),(N+2)*(Q+1));
Bm(1:Q+1,1:Q+1) = eye(Q+1);
Bm(Q+2:end,end-Q:end) = eye(Q+1);

Bf = [];
for i = 1:N+1
    Bf = blkdiag(Bf,[1 zeros(1,Q+1);zeros(1,Q+1) 1]);
end

B = blkdiag(Bm,Bf);

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

% matrices projection sur C %Cst
A = [D ; B]; 
delta = A*A'; 

normalise = @(f) f/sum(f(:));

sigma = 0.005; mini = 0.000001;
f0 = normalise(mini + gauss(0.1,sigma,N)); 
f1 = normalise(mini + gauss(0.9,sigma,N));

y = [zeros((N+1)*(Q+1),1) ; zeros(2*(Q+1),1) ; reshape([f1;f0],2*(N+1),1)];
Cst = A'*(delta\y);

P = eye((N+1)*(Q+2)+(N+2)*(Q+1)) - A'*(delta\A);
