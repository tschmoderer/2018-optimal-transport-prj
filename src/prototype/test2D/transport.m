clc
clear all
close all

globals; 

% %% Initialisation %%

N = 1;
P = 1;
Q = 1;

divergence;
boundary;
interpolation;

% matrices projection sur C %
 A = [D ; B]; 
 delta = A*A'; 
 
 normalise = @(f) f/sum(f(:));
 
 sigma = 0.05; mini = 0.0001;
 f0 = normalise(gauss(0.1,0.1,sigma,sigma,N,P,mini));
 f1 = normalise(gauss(0.1,0.1,sigma,sigma,N,P,mini));

% y = [zeros((N+1)*(Q+1),1) ; zeros(2*(Q+1),1) ; reshape([f1;f0],2*(N+1),1)];
% Cst = A'*(delta\y);
 
% P = eye((N+1)*(Q+2)+(N+2)*(Q+1)) - A'*(delta\A);

% 
% %% test %%
% 
% N = 1; P = 7; Q = 9;
% 
% m = zeros(N+1,P+1,Q+1,2); % c'est un champ de vecteurs sur le cube espace temps
% f = zeros(N+1,P+1,Q+1); % c'est un champ dans le cube
% 
% mbarX = zeros(N+2,P+1,Q+1,2); mbarY = zeros(N+1,P+2,Q+1,2);
% fbar  = zeros(N+1,P+1,Q+2);
% 
% % test interpolation 
% Im = 0.25*(mbarX(1:end-1,:,:,:) + mbarX(2:end,:,:,:) + mbarY(:,1:end-1,:,:) + mbarY(:,2:end,:,:));
% If = 0.5*(fbar(:,:,1:end-1) + fbar(:,:,2:end));
% 
% % test divergence
% d = N*(mbarX(2:end,:,:,1) - mbarX(1:end-1,:,:,1)) + P*(mbarY(:,2:end,:,2) - mbarY(:,1:end-1,:,2)) + Q*(fbar(:,:,2:end) - fbar(:,:,1:end-1));
% 
% % test boundary 
% bX1 = mbarX(1,:,:,:); bX2 = mbarX(end,:,:,:);
% bY1 = mbarY(:,1,:,:); bY2 = mbarY(:,end,:,:);
% bt1 = fbar(:,:,1); bt2 = fbar(:,:,end);
% 
% %% fin tests %%

m = zeros(N+1,P+1,Q+1,2);
f = zeros(N+1,P+1,Q+1);

V = [reshape(m,[],1); reshape(f,[],1)];

mbar = zeros(N+2,P+2,Q+1,2);
fbar = zeros(N+1,P+1,Q+2);

U = [reshape(mbar,[],1); reshape(fbar,[],1)];






















% %% Fin Initialisation %%
% 
% m = zeros(Q+1,N+1);
% f = zeros(Q+1,N+1);
% mbar = zeros(Q+1,N+2);
% fbar = zeros(Q+2,N+1);
% % %% Fin Initialisation %%
% 
% m = zeros(Q+1,N+1);
% f = zeros(Q+1,N+1);
% mbar = zeros(Q+1,N+2);
% fbar = zeros(Q+2,N+1);
% 
% V = [reshape(m,(N+1)*(Q+1),1);reshape(f,(N+1)*(Q+1),1)];
% U = [r% %% Initialisation %%
% 
% N = 32;
% P = 32;
% Q = 32;
% 
% % Matrice de l'opérateur b %
% Bm = zeros(2*(Q+1),(N+2)*(Q+1));
% Bm(1:Q+1,1:Q+1) = eye(Q+1);
% Bm(Q+2:end,end-Q:end) = eye(Q+1);
% 
% Bf = [];
% for i = 1:N+1
%     Bf = blkdiag(Bf,[1 zeros(1,Q+1);zeros(1,Q+1) 1]);
% end
% 
% B = blkdiag(Bm,Bf);
% % matrice d'interpolation %
% 
% Im = zeros((N+1)*(Q+1),(N+2)*(Q+1));
% for i = 1:(N+1)*(Q+1)
%     for j = 1:(N+2)*(Q+1)
%         if i == j 
%             Im(i,j) = 1;
%         elseif j == i+Q+1
%             Im(i,j) = 1;
%         end
%     end
% end
% dia = zeros(Q+1,Q+2);
% for i = 1:Q+1
%     for j = 1:Q+2
%         if i == j 
%             dia(i,j) = 1;
%         elseif j == i+1
%             dia(i,j) = 1;
%         end
%     end
% end
% If = [];
% for i = 1:N+1
%     = zeros(Q+1,Q+2);
% for i = 1:Q+1
%     for j = 1:Q+2
%         if i == j 
%             dia(i,j) = -1;
%         elseif j == i+1
%             dia(i,j) = 1;
%         end
%     end
% end
% Df = [];
% for i = 1:N+1
%     Df = blkdiag(Df,dia);
% end
% 
% D = [N*Dm Q*Df]; If = blkdiag(If,dia);
% end
% 
% Interp = 0.5*blkdiag(Im,If);
% % matrice de projection sur G2
% pG2 = inv(eye((N+1)*(Q+2)+(N+2)*(Q+1)) + Interp'*Interp);
% 
% % Matrice de la divergence %
% Dm = zeros((N+1)*(Q+1),(N+2)*(Q+1));
% for i = 1:(N+1)*(Q+1)
%     for j = 1:(N+2)*(Q+1)
%         if i == j 
%             Dm(i,j) = -1;
%         elseif j == i+Q+1
%             Dm(i,j) = 1;
%         end
%     end
% end
% dia = zeros(Q+1,Q+2);
% for i = 1:Q+1
%     for j = 1:Q+2
%         if i == j 
%             dia(i,j) = -1;
%         elseif j == i+1
%             dia(i,j) = 1;
%         end
%     end
% end
% Df = [];
% for i = 1:N+1
%     Df = blkdiag(Df,dia);
% end
% 
% D = [N*Dm Q*Df];
%reshape(mbar,(N+2)*(Q+1),1);reshape(fbar,(N+1)*(Q+2),1)];
% 
% sigma = 0.05; mini = 0.0001;
% f0 = gauss(0.2,sigma,N,mini); 
% f1 = gauss(0.8,sigma,N,mini); 
% 
% alpha = 0.5; gamma = 1;
% 
% wU0 = zeros(size(U)); wV0 = zeros(size(V));
% zU0 = zeros(size(U)); zV0 = zeros(size(V));
% 
% [XX,YY] = meshgrid([0:N],[0:Q]); YY = flipud(YY);
% 
% % Itérations
% niter = 500;
% cout = zeros(1,niter);
% div = zeros(1,niter);
% for l = 1:niter
%     [wU1 , wV1] = proxG1(2*zU0-wU0,2*zV0-wV0,gamma);
%     wU1 = wU0 + alpha*(wU1- zU0); wV1 = wV0 + alpha*(wV1- zV0);
%     
%     [zU0,zV0] = proxG2(wU1,wV1);
%     wU0 = wU1;
%     wV0 = wV1;
%     
%     f = reshape(zV0((N+1)*(Q+1)+1:end),Q+1,N+1);
%     surf(XX,YY,f)
%     xlabel('x');
%     ylabel('t');
%     zlabel('f');
%     title(['itération : ',num2str(l)]);
%     p% %% Initialisation %%
% 
% N = 32;
% P = 32;
% Q = 32;
% 
% % Matrice de l'opérateur b %
% Bm = zeros(2*(Q+1),(N+2)*(Q+1));
% Bm(1:Q+1,1:Q+1) = eye(Q+1);
% Bm(Q+2:end,end-Q:end) = eye(Q+1);
% 
% Bf = [];
% for i = 1:N+1
%     Bf = blkdiag(Bf,[1 zeros(1,Q+1);zeros(1,Q+1) 1]);
% end
% 
% B = blkdiag(Bm,Bf);
% % matrice d'interpolation %
% 
% Im = zeros((N+1)*(Q+1),(N+2)*(Q+1));
% for i = 1:(N+1)*(Q+1)
%     for j = 1:(N+2)*(Q+1)
%         if i == j 
%             Im(i,j) = 1;
%         elseif j == i+Q+1
%             Im(i,j) = 1;
%         end
%     end
% end
% dia = zeros(Q+1,Q+2);
% for i = 1:Q+1
%     for j = 1:Q+2
%         if i == j 
%             dia(i,j) = 1;
%         elseif j == i+1
%             dia(i,j) = 1;
%         end
%     end
% end
% If = [];
% for i = 1:N+1
%     = zeros(Q+1,Q+2);
% for i = 1:Q+1
%     for j = 1:Q+2
%         if i == j 
%             dia(i,j) = -1;
%         elseif j == i+1
%             dia(i,j) = 1;
%         end
%     end
% end
% Df = [];
% for i = 1:N+1
%     Df = blkdiag(Df,dia);
% end
% 
% D = [N*Dm Q*Df]; If = blkdiag(If,dia);
% end
% 
% Interp = 0.5*blkdiag(Im,If);
% % matrice de projection sur G2
% pG2 = inv(eye((N+1)*(Q+2)+(N+2)*(Q+1)) + Interp'*Interp);
% 
% % Matrice de la divergence %
% Dm = zeros((N+1)*(Q+1),(N+2)*(Q+1));
% for i = 1:(N+1)*(Q+1)
%     for j = 1:(N+2)*(Q+1)
%         if i == j 
%             Dm(i,j) = -1;
%         elseif j == i+Q+1
%             Dm(i,j) = 1;
%         end
%     end
% end
% dia = zeros(Q+1,Q+2);
% for i = 1:Q+1
%     for j = 1:Q+2
%         if i == j 
%             dia(i,j) = -1;
%         elseif j == i+1
%             dia(i,j) = 1;
%         end
%     end
% end
% Df = [];
% for i = 1:N+1
%     Df = blkdiag(Df,dia);
% end
% 
% D = [N*Dm Q*Df];
%pause(0.04)
%     
%     cout(i) = cost(zV0);
%     div(i) = sum(D*zU0);
% end
% 
% figure;
% subplot(2,1,1)
% plot([1:niter],cout);
% title('cout')
% subplot(2,1,2)
% plot([1:niter],div)
% title('div')
% V = [reshape(m,(N+1)*(Q+1),1);reshape(f,(N+1)*(Q+1),1)];
% U = [reshape(mbar,(N+2)*(Q+1),1);reshape(fbar,(N+1)*(Q+2),1)];
% 
% sigma = 0.05; mini = 0.0001;
% f0 = gauss(0.2,sigma,N,mini); 
% f1 = gauss(0.8,sigma,N,mini); 
% 
% alpha = 0.5; gamma = 1;
% 
% wU0 = zeros(size(U)); wV0 = zeros(size(V));
% zU0 = zeros(size(U)); zV0 = zeros(size(V));
% 
% [XX,YY] = meshgrid([0:N],[0:Q]); YY = flipud(YY);
% 
% % Itérations
% niter = 500;
% cout = zeros(1,niter);
% div = zeros(1,niter);
% for l = 1:niter
%     [wU1 , wV1] = proxG1(2*zU0-wU0,2*zV0-wV0,gamma);
%     wU1 = wU0 + alpha*(wU1- zU0); wV1 = wV0 + alpha*(wV1- zV0);
%     
%     [zU0,zV0] = proxG2(wU1,wV1);
%     wU0 = wU1;
%     wV0 = wV1;
%     
%     f = reshape(zV0((N+1)*(Q+1)+1:end),Q+1,N+1);
%     surf(XX,YY,f)
%     xlabel('x');
%     ylabel('t');
%     zlabel('f');
%     title(['itération : ',num2str(l)]);
%     pause(0.04)
% %% Initialisation %%
% 
% N = 32;
% P = 32;
% Q = 32;
% 
% % Matrice de l'opérateur b %
% Bm = zeros(2*(Q+1),(N+2)*(Q+1));
% Bm(1:Q+1,1:Q+1) = eye(Q+1);
% Bm(Q+2:end,end-Q:end) = eye(Q+1);
% 
% Bf = [];
% for i = 1:N+1
%     Bf = blkdiag(Bf,[1 zeros(1,Q+1);zeros(1,Q+1) 1]);
% end
% 
% B = blkdiag(Bm,Bf);
% % matrice d'interpolation %
% 
% Im = zeros((N+1)*(Q+1),(N+2)*(Q+1));
% for i = 1:(N+1)*(Q+1)
%     for j = 1:(N+2)*(Q+1)
%         if i == j 
%             Im(i,j) = 1;
%         elseif j == i+Q+1
%             Im(i,j) = 1;
%         end
%     end
% end
% dia = zeros(Q+1,Q+2);
% for i = 1:Q+1
%     for j = 1:Q+2
%         if i == j 
%             dia(i,j) = 1;
%         elseif j == i+1
%             dia(i,j) = 1;
%         end
%     end
% end
% If = [];
% for i = 1:N+1
%     = zeros(Q+1,Q+2);
% for i = 1:Q+1
%     for j = 1:Q+2
%         if i == j 
%             dia(i,j) = -1;
%         elseif j == i+1
%             dia(i,j) = 1;
%         end
%     end
% end
% Df = [];
% for i = 1:N+1
%     Df = blkdiag(Df,dia);
% end
% 
% D = [N*Dm Q*Df]; If = blkdiag(If,dia);
% end
% 
% Interp = 0.5*blkdiag(Im,If);
% % matrice de projection sur G2
% pG2 = inv(eye((N+1)*(Q+2)+(N+2)*(Q+1)) + Interp'*Interp);
% 
% % Matrice de la divergence %
% Dm = zeros((N+1)*(Q+1),(N+2)*(Q+1));
% for i = 1:(N+1)*(Q+1)
%     for j = 1:(N+2)*(Q+1)
%         if i == j 
%             Dm(i,j) = -1;
%         elseif j == i+Q+1
%             Dm(i,j) = 1;
%         end
%     end
% end
% dia = zeros(Q+1,Q+2);
% for i = 1:Q+1
%     for j = 1:Q+2
%         if i == j 
%             dia(i,j) = -1;
%         elseif j == i+1
%             dia(i,j) = 1;
%         end
%     end
% end
% Df = [];
% for i = 1:N+1
%     Df = blkdiag(Df,dia);
% end
% 
% D = [N*Dm Q*Df];
%     
%     cout(i) = cost(zV0);
%     div(i) = sum(D*zU0);
% end
% 
% figure;
% subplot(2,1,1)
% plot([1:niter],cout);
% title('cout')
% subplot(2,1,2)
% plot([1:niter],div)
% title('div')