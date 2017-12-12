clc
clear all
close all

%% test poisson -- OK
% 
% % N = 80; Q = 50; 
% % X = [0:N]/N; Y = [0:Q]/Q;
% % [XX,YY] = meshgrid(X,Y); %YY = flipud(YY);
% % 
% % bL = Y'.^2; bR = 1+Y'.^2;
% % bU = 1+X.^2; bD = X.^2;
% % d = -4*ones(Q+1,N+1);
% % %bL = 0; bR = 0; bU = 0; bD = 1;
% % S = poisson(d,bL,bR,bU,bD,N,Q,1e-15);
% % %S = poisson(-ones(Q+1,N+1),0,0,0,0,N,Q,1e-5);
% % % X = [0:N]/N; Y = [0:Q]/Q;
% % %S = poisson(-6*XX.*YY.*(1-YY)+2*XX.^3,0,Y.*(1-Y),0,0,N,Q,1e-5);
% % 
% % surf(XX,YY,S)
% % %sum(sum(S-YY.*(1-YY).*XX.^3))
% % sum(sum(abs(XX.^2+YY.^2 - S)))
% % 
% % errL = sum(abs(S(:,1) - bL));
% % errR = sum(abs(S(:,end) - bR));
% % errU = sum(abs(S(1,:) - bU));
% % errD = sum(abs(S(end,:) - bD));
% % [errL errR errU errD]
% % 
% % return


N = 11;
Q = 9; 

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

D = [N*Dm Q*Df];%!!

% matrices projection sur C %
A = [D ; B]; 
delta = A*A'; 


normalise = @(f) f/sum(f(:));

f0 = normalise(gauss(0.4,0.005,N) + 0.01);
f1 = normalise(gauss(0.6,0.005,N) + 0.01);

yy = zeros(Q+1,N+1);
yy(1,:) = Q*f1; 
yy(end,:) = -Q*f0;
%C1 = poisson(yy,0,0,0,0,N,Q,1e-10);
%C1 = poisson(zeros(Q+1,N+1),0,0,-Q*f1,Q*f0,N,Q,1e-10);
C1 = poisson2d_Neumann(yy);
Cmbar = zeros(Q+1,N+2); Cfbar = zeros(Q+2,N+1);
    
Cmbar(:,1) = -C1(:,1); Cmbar(:,end) = C1(:,end);  Cmbar(:,1) = 0; Cmbar(:,end) = 0;
Cmbar(:,2:end-1) = (C1(:,1:end-1) - C1(:,2:end));

Cfbar(1,:) = -C1(1,:); Cfbar(end,:) = C1(end,:); Cfbar(1,:) = f1/Q; Cfbar(end,:) = f0/Q;
Cfbar(2:end-1,:) = (C1(1:end-1,:) - C1(2:end,:));

Cmbar = N*Cmbar; Cfbar = Q*Cfbar;

y = [zeros((N+1)*(Q+1),1) ; zeros(2*(Q+1),1) ; reshape([f1;f0],2*(N+1),1)];
C2 = delta\y;
C22 = A'*C2;

errCs = sum(sum(abs(reshape(C2(1:(N+1)*(Q+1)),Q+1,N+1) - C1)))
errCm = sum(sum(abs(reshape(C22(1:(N+2)*(Q+1)),Q+1,N+2) - Cmbar)))
eerCf = sum(sum(abs(reshape(C22((N+2)*(Q+1)+1:end),Q+2,N+1) - Cfbar)))


figure;

subplot(2,2,1), 
surf(Cmbar);
title('Cmbar');

subplot(2,2,2), 
surf(Cfbar);
title('Cfbar');

subplot(2,2,3), 
surf(reshape(C22(1:(N+2)*(Q+1)),Q+1,N+2));
title('C2 - mbar');

subplot(2,2,4), 
surf(reshape(C22((N+2)*(Q+1)+1:end),Q+2,N+1));
title('C2 - fbar');

figure; 
subplot(121)
surf(C1);
title('C1 - poisson');
subplot(122);
surf(reshape(C2(1:(N+1)*(Q+1)),Q+1,N+1))
title('C2 - delta\y');

sum(sum(abs(C1-reshape(C2(1:(N+1)*(Q+1)),Q+1,N+1))))

return
% f = zeros(Q+1,N+1); f(1,:) = -Q*f1; f(end,:) = Q*f0; 
% surf(poisson2d_Neumann(f))


























%% test opérateur 
mbar = rand(Q+1,N+2);
fbar = rand(Q+2,N+1); 

U = [reshape(mbar,(N+2)*(Q+1),1);reshape(fbar,(N+1)*(Q+2),1)];

d = N*(mbar(:,2:end)-mbar(:,1:end-1))+Q*(fbar(2:end,:)-fbar(1:end-1,:));

dU = D*U;
sum(sum(abs(d-reshape(dU,Q+1,N+1))))

bU = B*U;
b = [mbar(:,1); mbar(:,end) ; reshape([fbar(1,:);fbar(end,:)],2*(N+1),1)];
sum(sum(abs(bU-b)))


C1 = rand(Q+1,N+1); 
Cmbar = zeros(Q+1,N+2); Cfbar = zeros(Q+2,N+1);
    
Cmbar(:,1) = -C1(:,1); Cmbar(:,end) = C1(:,end);
Cmbar(:,2:end-1) = C1(:,1:end-1) - C1(:,2:end);

Cfbar(1,:) = -C1(1,:); Cfbar(end,:) = C1(end,:); 
Cfbar(2:end-1,:) = C1(1:end-1,:) - C1(2:end,:);

Cmbar = N*Cmbar; Cfbar = Q*Cfbar;

dpV = D'*reshape(C1,[],1);



sum(sum(abs(Cmbar - reshape(dpV(1:(N+2)*(Q+1)),Q+1,N+2))))
sum(sum(abs(Cfbar - reshape(dpV((N+2)*(Q+1)+1:end),Q+2,N+1))))