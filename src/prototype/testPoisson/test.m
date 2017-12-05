clc
clear all
close all

N = 31;
Q = 29; 

% Matrice de l'opérateur b
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

%D = [N*Dm Q*Df];%!!
D = [Dm Df];
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
% 
figure;
subplot(2,2,1), 
surf(Cmbar);
%contour3(C1,50);
title('Cmbar');
subplot(2,2,2), 
surf(Cfbar);
%contour3(C1,50);
title('Cfbar');
subplot(2,2,3), 
surf(reshape(C22(1:(N+2)*(Q+1)),Q+1,N+2));
%contour3(reshape(C2(1:(N+1)*(Q+1)),Q+1,N+1),50);
title('C2 - mbar');
subplot(2,2,4), 
surf(reshape(C22((N+2)*(Q+1)+1:end),Q+2,N+1));
%contour3(reshape(C2(1:(N+1)*(Q+1)),Q+1,N+1),50);
title('C2 - fbar');

figure; 
subplot(121)
surf(C1);
title('C1 - poisson');
subplot(122);
surf(reshape(C2(1:(N+1)*(Q+1)),Q+1,N+1))
title('C2 - delta\y');
%close all 

% N = 10; P = 10;
% w = ([0:10]'*ones(1,P+1))'
% M = zeros(P+1);



% %% test div.div*
% N = 10; Q = 10; 
% 
% w = rand(Q+1,N+1); % en theorie les inconnue
% 
% Lw = zeros(Q+1,N+1);
% 
% Lw(2:end-1,2:end-1) = (-w(2:end-1,1:end-2)+2*w(2:end-1,2:end-1)-w(2:end-1,3:end)) + (-w(1:end-2,2:end-1)+2*w(2:end-1,2:end-1)-w(3:end,2:end-1));
% 
% Lw(2:end-1,1) = 2*w(2:end-1,1) - w(2:end-1,2) - w(1:end-2,1) + 2*w(2:end-1,1) - w(3:end,1); % bord gauche
% Lw(2:end-1,end) = 2*w(2:end-1,end) - w(2:end-1,end-1) - w(1:end-2,end) + 2*w(2:end-1,end) - w(3:end,end); % bord droit
% 
% Lw;







%% test subroutine

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