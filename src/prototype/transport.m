clc
close all
clear

%%% Initialisation %%%
globals;

d = 1; % dimension du système

N = 200; % Nb de points de discrétisation dans le sens des x
%P = 5; % Nb de points de discrétisation dans le sens des y
Q = 100; % Nb de points de discrétisation dans le sens de t

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

m = rand(size(m));
f = rand(size(f));
mbar = rand(size(mbar));
fbar = rand(size(fbar));

f0 = f0(Gc.x);
f1 = f1(Gc.x);


% Ecrire les valeurs dans les fichiers.
fiD = fopen('files/parameters.txt','w');
f0Id = fopen('files/f0.txt','w');
f1Id = fopen('files/f1.txt','w');

fprintf(fiD,'%d %d',N,Q); % FreeFem needs 
fprintf(f0Id,'%12f\n',f0);
fprintf(f1Id,'%12s\r\n',f1);

fclose(fiD);
fclose(f0Id);
fclose(f1Id);


%%% Zone de tests %%%

% test calcul de la constante %%
y = dlmread('files/Y.txt'); % size NxQ
y = reshape(y,N+1,Q+1)';

[X,Y] = meshgrid(Gc.x,Gc.t);
surf(X,Y,y);
xlabel('x')
ylabel('t');

%%% Fin zone de tests %%%