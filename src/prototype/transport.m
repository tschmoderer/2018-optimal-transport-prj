clc
close all
clear

%%% Initialisation %%%
globals;

N = 29; % Nb de points de discrétisation dans le sens des x
Q = 33; % Nb de points de discrétisation dans le sens de t

%% Grille centrée %%
Gc.x = [0:N]/N;
Gc.t = [0:Q]/Q;

%% Grille décentrée %%
Gs.x = ([-1:N]+0.5)/N;
Gs.t = ([-1:Q]+0.5)/Q;

% la grille de temps : 
t = repmat(linspace(1,0,Q+1)',1,N+1);

%% Variables centrées %%
m = zeros(Q+1,N+1);
f = zeros(Q+1,N+1);

%% Variables décentrées %%
mbar = zeros(Q+1,N+2);
fbar = zeros(Q+2,N+1);

%% Frontières %%
b0.m = [zeros(Q+1,1) , zeros(Q+1,1)];
b0.f = [f1(Gc.x) ; f0(Gc.x)];

%% Matrice d'interpolation %%
Interpm = [diag(ones(1,N+1));zeros(1,N+1)] + [zeros(1,N+1);diag(ones(1,N+1))];
Interpm = Interpm/2; % m = mbar*I
Interpf = [diag(ones(1,Q+1)) zeros(Q+1,1)] + [zeros(Q+1,1) diag(ones(1,Q+1))];
Interpf = Interpf/2;

Interpm_adj = Interpm';
Interpf_adj = Interpf';

%% Paramètres %%
alpha = 1.983; % Doit etre dans ]0,2[
beta = 1; % Doit etre dans [0,1]
gamma = 1/240; % Doit etre positif

%%% Fin Initialisation %%%

[XX YY] = meshgrid(linspace(0,1,N+1),linspace(0,1,Q+1));
YY  = flipud(YY);
[XXX YYY] = meshgrid(linspace(0,1,N+2),linspace(0,1,Q+1));
[XXXX YYYY] = meshgrid(linspace(0,1,N+1),linspace(0,1,Q+2));

%system("FreeFem++ -v 0 poisson_2d_constante.pde");
Cst = poisson(zeros(Q+1,N+1),0,0,b0.f(1,:),b0.f(2,:),N,Q,1e-3);
[Cstmbar, Cstfbar] = divergence_adjoint(Cst);

% % Initialisation 
 Zm = zeros(Q+1,N+1); Zf = zeros(Q+1,N+1);
 Wm0 = zeros(Q+1,N+1); Wf0 = zeros(Q+1,N+1); Wm1 = zeros(Q+1,N+1); Wf1 = zeros(Q+1,N+1); Wm2 = zeros(Q+1,N+1); Wf2 = zeros(Q+1,N+1);
 Zmbar = zeros(Q+1,N+2); Zfbar = zeros(Q+2,N+1); 
 Wmbar0 = zeros(Q+1,N+2); Wfbar0 = zeros(Q+2,N+1); Wmbar1 = zeros(Q+1,N+2); Wfbar1 = zeros(Q+2,N+1); Wmbar2 = zeros(Q+1,N+2); Wfbar2 = zeros(Q+2,N+1); 

 Wf0 = (1-t).*repmat(b0.f(2,:),Q+1,1) + t.*repmat(b0.f(1,:),Q+1,1);
% % Itérations
 niter = 250;
 cout = zeros(1,niter);
for i = 1:niter
	[Zmbar,Zfbar,Zm,Zf] = proxG2(Wmbar0,Wfbar0,Wm0,Wf0);

 	[Wmbar1,Wfbar1,Wm1,Wf1] = proxG1(Wmbar0,Wfbar0,Wm0,Wf0,gamma);
     
 	Wmbar1 = 2*Wmbar1 - Wmbar0; Wfbar1 = 2*Wfbar1 - Wfbar0;
 	Wm1 = 2*Wm1 - Wm0; Wf1 = 2*Wf1 - Wf0;	
 
 	[Wmbar2,Wfbar2,Wm2,Wf2] = proxG2(Wmbar1,Wfbar1,Wm1,Wf1);
 
    Wmbar2 = 2*Wmbar2 - Wmbar1; Wfbar2 = 2*Wfbar2-Wfbar1; 
 	Wm2 = 2*Wm2-Wm1; Wf2=2*Wf2-Wf1;
 	
 	Wmbar0 = (1-0.5*alpha)*Wmbar0 + 0.5*alpha*Wmbar2;
 	Wfbar0 = (1-0.5*alpha)*Wfbar0 + 0.5*alpha*Wfbar2;
 	Wm0 = (1-0.5*alpha)*Wm0 + 0.5*alpha*Wm2;
 	Wf0 = (1-0.5*alpha)*Wf0 + 0.5*alpha*Wf2;
     
	surf(XX,YY,Zf)
	xlabel('x')
	ylabel('t')
    title(['iteration : ',num2str(i)])
	drawnow
%    pause
    
    cout(i) = cost(Zm,Zf);
end

plot([1:niter],cout);
title('cout');
xlabel('iter');
ylabel('cout');
