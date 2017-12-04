clc
clear all
close all

N = 29;
Q = 31; 

% Matrice de l'opérateur b
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

% matrice de projection sur G2 %
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


normalise = @(f) f/sum(f(:));

f0 = normalise(gauss(0.4,0.05,N) + 0.01);
f1 = normalise(gauss(0.6,0.05,N) + 0.01);

C1 = poisson(zeros(Q+1,N+1),0,0,f1,f0,N,Q,1e-10);
Cmbar = zeros(Q+1,N+2); Cfbar = zeros(Q+2,N+1);
    
Cmbar(:,1) = -C1(:,1); Cmbar(:,end) = C1(:,end);
Cmbar(:,2:end-1) = C1(:,1:end-1) - C1(:,2:end);

Cfbar(1,:) = -C1(1,:); Cfbar(end,:) = C1(end,:); 
Cfbar(2:end-1,:) = C1(1:end-1,:) - C1(2:end,:);

Cmbar = Cmbar; Cfbar = Cfbar;

y = [zeros((N+1)*(Q+1),1) ; zeros(2*(Q+1),1) ; reshape([f1;f0],2*(N+1),1)];
C2 = delta\y;
C22 = A'*C2;

errCm = sum(sum(reshape(C22(1:(N+2)*(Q+1)),Q+1,N+2) - Cmbar))
eerCf = sum(sum(reshape(C22((N+2)*(Q+1)+1:end),Q+2,N+1) - Cfbar))

figure;
subplot(1,2,1), 
surf(C1);
%contour3(C1,50);
title('C1');

subplot(1,2,2), 
surf(reshape(C2(1:(N+1)*(Q+1)),Q+1,N+1));
%contour3(reshape(C2(1:(N+1)*(Q+1)),Q+1,N+1),50);
title('C2');


% m = rand(Q+1,N+1); 
% f = rand(Q+1,N+1);
% 
% mbar = rand(Q+1,N+2);
% fbar = rand(Q+2,N+1);
% 
% cost(m,f)
% proxJ(m,f)
% 
% proxG2(mbar,fbar,m,f,N,Q)
% 
% projC(mbar,fbar,N,Q)
% 
% surf(ans)