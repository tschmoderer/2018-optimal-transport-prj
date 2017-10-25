clc
close all
clear

%%% Initialisation %%%

d = 1; % dimension du système

N = 5; % Nb de points de discrétisation dans le sens des x
%P = 5; % Nb de points de discrétisation dans le sens des y
Q = 2; % Nb de points de discrétisation dans le sens de t

%% Grille centrée %%
Gc.x = [0:N]/N;
%Gc.y = [0:P]/P;
Gc.t = [0:Q]/Q;

%% Grille décentrée %%
Gs.x = ([-1:N]+0.5)/N;
%Gs.y = ([-1:P+1]+0.5)/P;
Gs.t = ([-1:Q]+0.5)/Q;

%% Variables centrées %%
m = zeros(Q+1,N+1);
f = zeros(Q+1,N+1);

%% Variables décentrées %%
mbar = zeros(Q+1,N+2);
fbar = zeros(Q+2,N+1);

%% Frontières %%
b0.m = [zeros(Q+1,1) , zeros(Q+1,1)];
b0.f = [f0(Gc.x) ; f1(Gc.x)];

%% Matrice d'interpolation %%
Interpm = [diag(ones(1,N+1));zeros(1,N+1)] + [zeros(1,N+1);diag(ones(1,N+1))];
Interpm = Interpm/2; % m = mbar*I

Interpm_adj = Interpm';

Interpf = [diag(ones(1,Q+1)) zeros(Q+1,1)] + [zeros(Q+1,1) diag(ones(1,Q+1))];
Interpf = Interpf/2;

Interpf_adj = Interpf';

%% Matrice de divergence %%

%% Matrice d'extraction des frontières %%


%% Paramètres %%
alpha = 1; % Doit etre dans ]0,2[
beta = 1; % Doit etre dans [0,1]
gamma = 0.5; % Doit etre positif

%%% Fin Initialisation %%%


%%% Zone de tests %%%


%%% Fin zone de tests %%%
m = rand(size(m));
f = rand(size(f));
mbar = rand(size(mbar));
fbar = rand(size(fbar));


f0 = f0(Gc.x);
f1 = f1(Gc.x);

fiD = fopen('parameters.txt','w');
f0Id = fopen('f0.txt','w');
f1Id = fopen('f1.txt','w');

fprintf(fiD,'%d\r\n%d',N,Q);
fprintf(f0Id,'%12s\r\n',f0);
fprintf(f1Id,'%12s\r\n',f1);

fclose(fiD);
fclose(f0Id);
fclose(f1Id);

dlmread("Y.txt");

[X,Y] = meshgrid(linspace(0,1,501),linspace(0,1,101));
test = dlmread("test.txt"); % size 21*11
test = reshape(test,101,501);
surf(X,Y,test);