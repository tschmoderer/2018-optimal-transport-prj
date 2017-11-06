% variables globales 
global N; % Nb de points de discrétisation dans le sens des x
%P = 5; % Nb de points de discrétisation dans le sens des y
global Q; % Nb de points de discrétisation dans le sens de t

% %% Grille centrée %%
% global Gc;
% 
% %% Grille décentrée %%
% global Gs;

%%% Variables centrées %%
%global m;
%global f;
%
%%% Variables décentrées %%
%global mbar;
%global fbar;

%% Frontières %%
global b0;

%% Matrice d'interpolation %%
global Interpm;
global Interpm_adj;

global Interpf;
global Interpf_adj;

%% Paramètres %%
global alpha; % Doit etre dans ]0,2[
global beta; % Doit etre dans [0,1]
global gamma; % Doit etre positif

%% Constantes dans la projection sur C %%
global Cst;
global Cstmbar;
global Cstfbar;