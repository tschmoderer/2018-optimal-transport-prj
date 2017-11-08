clc
clear
close all 

% parametres
N = 20; P = 20;
alpha = 1.983; gamma = 1/240;
niter = 50;

% Données 
X = linspace(0,1,N+1);
minimal = 0.05;
gauss = @(mu,sigma)(minimal + exp(-0.5*((X-mu)/sigma).^2));
normal = @(f)(f/sum(f(:)));

f0 = normal(gauss(0.1,0.05))';
f1 = normal(gauss(0.9,0.05))';

b0.m = [zeros(1,P+1);zeros(1,P+1)];
b0.f = [f0,f1];

cost = @(w)(sum(w.m(:).^2./w.f(:)));

% Polynome dans le proximal de J 

















t = repmat(linspace(0,1,P+1),N+1,1);

w0.m = zeros(N+1,P+1);
w0.f = (1-t).*repmat(f0,1,P+1) + t.*repmat(f1,1,P+1);

[XX YY] = meshgrid(linspace(0,1,N+1),linspace(0,1,P+1));
surf(XX,YY,w0.f);
xlabel('t'); ylabel('x');
pause

RproxG1 = @(w)(2*ProxG1(w) - w);
RproxG2 = @(w)(2*ProxG2(w) - w);

% % Boucle principale 
% for i = 1:niter
%    w1 = (1-0.5*alpha)*w0 + 0.5*alpha*RproxG2(RproxG1(w0)); 
% end