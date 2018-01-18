clc
clear all
close all

N = 30;
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

D = [N*Dm Q*Df];

% matrices projection sur C %
A = [D ; B]; 
delta = A*A'; 

normalise = @(f) f/sum(f(:));

f0 = normalise(gauss(0.4,0.005,N) + 0.01);
f1 = normalise(gauss(0.6,0.005,N) + 0.01);

%!!%
ymbar = zeros(Q+1,N+2); 
yfbar = zeros(Q+2,N+1); yfbar(1,:) = f1; yfbar(end,:) = f0;
yy = N*(ymbar(:,2:end)-ymbar(:,1:end-1)) + Q*(yfbar(2:end,:) - yfbar(1:end-1,:));

C1 = poisson2d_Neumann(-yy);
% C1 = poisson_neumann(zeros(Q+1,N+1),zeros(Q+1,1),zeros(Q+1,1),f1,f0,N,Q,1e-5);
% C1 = poisson_neumann(yy,zeros(Q+1,1),zeros(Q+1,1),0,0,N,Q,1e-5);

Cmbar = zeros(Q+1,N+2); Cfbar = zeros(Q+2,N+1);
    
Cmbar(:,1) = -C1(:,1); Cmbar(:,end) = C1(:,end);  Cmbar(:,1) = 0; Cmbar(:,end) = 0;
Cmbar(:,2:end-1) = (C1(:,1:end-1) - C1(:,2:end));

Cfbar(1,:) = -C1(1,:); Cfbar(end,:) = C1(end,:); Cfbar(1,:) = Cfbar(1,:) + f1/Q; Cfbar(end,:) = Cfbar(end,:) + f0/Q;
Cfbar(2:end-1,:) = (C1(1:end-1,:) - C1(2:end,:));

Cmbar = N*Cmbar; Cfbar = Q*Cfbar; 

y = [zeros((N+1)*(Q+1),1) ; zeros(2*(Q+1),1) ; reshape([f1;f0],2*(N+1),1)];
C2 = delta\y;
C22 = A'*C2;

errCs = sum(sum(abs(reshape(C2(1:(N+1)*(Q+1)),Q+1,N+1) - C1)))
errCm = sum(sum(abs(reshape(C22(1:(N+2)*(Q+1)),Q+1,N+2) - Cmbar)))
errCf = sum(sum(abs(reshape(C22((N+2)*(Q+1)+1:end),Q+2,N+1) - Cfbar)))


% !!%% 
% PBM 

mustBeCloseToZeros = sum(sum(abs(N*(Cmbar(:,2:end)-Cmbar(:,1:end-1)) + Q*(Cfbar(2:end,:) - Cfbar(1:end-1,:)))))


A = @(mbar,fbar) N*(Cmbar(:,2:end)-Cmbar(:,1:end-1)) + Q*(Cfbar(2:end,:) - Cfbar(1:end-1,:));


% 
% est non nulle
% !!%

% figure;
% 
% subplot(2,2,1), 
% surf(Cmbar);
% title('Cmbar');
% 
% subplot(2,2,2), 
% surf(Cfbar);
% title('Cfbar');
% 
% subplot(2,2,3), 
% surf(reshape(C22(1:(N+2)*(Q+1)),Q+1,N+2));
% title('C2 - mbar');
% 
% subplot(2,2,4), 
% surf(reshape(C22((N+2)*(Q+1)+1:end),Q+2,N+1));
% title('C2 - fbar');

figure; 
subplot(121),
surf(Cmbar);
hold on 
surf(reshape(C22(1:(N+2)*(Q+1)),Q+1,N+2),'EdgeColor','none');
title('Cmbar vs delta mbar (fade)');
hold off

subplot(122),
surf(Cfbar);
hold on 
surf(reshape(C22((N+2)*(Q+1)+1:end),Q+2,N+1),'EdgeColor','none');
title('Cfbar vs delta fbar (fade)');
hold off

% figure; 
% subplot(121)
% surf(C1);
% title('C1 - poisson');
% subplot(122);
% surf(reshape(C2(1:(N+1)*(Q+1)),Q+1,N+1))
% title('C2 - delta y');

figure; 
surf(C1);
hold on; 
surf(reshape(C2(1:(N+1)*(Q+1)),Q+1,N+1),'EdgeColor','none')
title('C1 vs Delta (fade)');

%DeltaSurC1 = reshape(C2(1:(N+1)*(Q+1)),Q+1,N+1)./C1


return



























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