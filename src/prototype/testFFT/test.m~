clc
%clear all
close all 

N = 2;
Q = 1;

m = rand(Q+1,N+1);
f = rand(Q+1,N+1); 

mbar = rand(Q+1,N+2);
fbar = rand(Q+2,N+1);

%% Test pour faire de la divergence un opérateur de convolution %%

dm = mbar(:,2:end) - mbar(:,1:end-1);
df = fbar(2:end,:) - fbar(1:end-1,:);

% 2D 
X = mbar;
Y = -[-1 1];
X1 = [X zeros(size(X,1),size(Y,2)-1);zeros(size(Y,1)-1,size(X,2)+size(Y,2)-1)];
Y1 = zeros(size(X1));    Y1(1:size(Y,1),1:size(Y,2)) = Y;
c1 = ifft2(fft2(X1).*fft2(Y1));
tm = conv2(X,Y,'full'); 

% 2D 
X = fbar;
Y = -[-1 1]';
X1 = [X zeros(size(X,1),size(Y,2)-1);zeros(size(Y,1)-1,size(X,2)+size(Y,2)-1)];
Y1 = zeros(size(X1)); Y1(1:size(Y,1),1:size(Y,2)) = Y;
c1 = ifft2(fft2(X1).*fft2(Y1));
tf = conv2(X,Y,'full'); 

eM = sum(sum(abs(dm-tm(:,2:end-1))))
eF = sum(sum(abs(df-tf(2:end-1,:))))


% pour evoiter l'extraction finale , je peux essayer de générer un filtre
% de taille (Q+1)*(N+1)

Dm = zeros(Q+1,N+1);
Dm(1,1) = 1; Dm(1,2) = -1;

conv2(Dm,mbar,'same')
dm
