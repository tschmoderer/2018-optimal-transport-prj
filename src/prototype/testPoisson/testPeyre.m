mynorm = @(a)norm(a(:));
sum3 = @(a)sum(a(:));

n = 20;
p = 20;

[Y,X] = meshgrid(linspace(0,1,n), linspace(0,1,n));
gaussian = @(a,b,sigma)exp( -((X-a).^2+(Y-b).^2)/(2*sigma^2) );
normalize = @(u)u/sum(u(:));

sigma = .1;
rho = .05; % minimum density value
f0 = normalize( rho + gaussian(.2,.3,sigma) );
f1 = normalize( rho + gaussian(.6,.7,sigma*.7) + .6*gaussian(.7,.4,sigma*.7) );

%%
% Display \(f_0,f_1\).

clf;
imageplot({f0 f1});

%%
% Boundary conditions, either periodic or Neumann.

bound = 'per';
%bound = 'neum';

%% 
% We use first order finite differences with periodic boundary condition to
% approximate the spacial derivatives.

if strcmp(bound, 'per')
    dx = @(u)u([2:end 1],:,:)-u;
    dy = @(u)u(:,[2:end 1],:)-u;
else
    dx = @(u)u([2:end end],:,:)-u;
    dy = @(u)u(:,[2:end end],:)-u;
end

%%
% The adjoint operators are backward derivatives.

if strcmp(bound, 'per')
    dxS = @(u)-u+u([end 1:end-1],:,:);
    dyS = @(u)-u+u(:,[end 1:end-1],:);
else
    dxS = @(u)[-u(1,:,:); u(1:end-2,:,:)-u(2:end-1,:,:); u(end-1,:,:)];
    dyS = @(u)[-u(:,1,:), u(:,1:end-2,:)-u(:,2:end-1,:), u(:,end-1,:)];    
end

%%
% Check that |dxS| and |dyS| really implement \(d/dx^*\) and \(d/dy^*\).

fprintf('Should be 0: %.2e\n', certify_adjoint(dx,dxS,[n n p]));
fprintf('Should be 0: %.2e\n', certify_adjoint(dy,dyS,[n n p]));

%%
% Define spacial gradient and divergence, satisfying \(\text{div}=-\nabla^*\).

grad = @(f)cat(4, dx(f), dy(f));
div  = @(u)-dxS(u(:,:,:,1)) - dyS(u(:,:,:,2));

%%
% Check that |div| really implements div.

fprintf('Should be 0: %.2e\n', certify_adjoint(grad,@(v)-div(v),[n n p]));

%%
% We use first order finite differences for the time derivatives. 
% Note that zero is padded at the end to keep the same dimensionality.

dt  = @(f)cat(3, f(:,:,2:end)-f(:,:,1:end-1), zeros(size(f,1),size(f,2)) );
dtS = @(u)cat(3, -u(:,:,1), u(:,:,1:end-2)-u(:,:,2:end-1), u(:,:,end-1));

A = @(w)cat( 3, div(w(:,:,:,1:2))+dt(w(:,:,:,3)), w(:,:,1,3), w(:,:,end,3) );

U = @(r0,r1)cat(3, r0, zeros(n,n,p-2), r1);
AS = @(s)cat(4, -grad(s(:,:,1:p)), dtS(s(:,:,1:p)) + U(s(:,:,end-1),s(:,:,end)) );

fprintf('Should be 0: %.2e\n', certify_adjoint(A,AS,[n n p 3]));

r0 = cat(3, zeros(n,n,p), f0, f1);

J = @(w)sum3(  sum(w(:,:,:,1:2).^2,4) ./ w(:,:,:,3)   );  

PolyCoef = @(m0,f0,lambda)[ones(length(f0),1), 4*lambda-f0, 4*lambda^2-4*f0, -lambda*sum(m0.^2,2) - 4*lambda^2*f0];

extract = @(A)A(:,1);
CubicReal = @(P)real( extract(poly_root(P')') );

Proxj0 = @(m0,f, lambda)cat(2, m0 ./ repmat( 1+2*lambda./f, [1 2]), f );
Proxj  = @(m0,f0,lambda)Proxj0( m0, CubicReal(PolyCoef(m0,f0,lambda)), lambda );

ProxJ = @(w,lambda)reshape( Proxj( ...
                   reshape(w(:,:,:,1:2), [n*n*p 2]), ...
                   reshape(w(:,:,:,3  ), [n*n*p 1]), lambda ), [n n p 3] ); 

%% Orthogonal Projection on the Constraints
% The proximal operator of the indicator function \(G=\iota_\Cc\) of \(\Cc\) is the
% projector, and does not depend on \(\la\).
% \[ \text{Prox}_{\gamma \iota_\Cc}(x)_i = 
%       \text{Proj}_\Cc(w) = w + A^* (A A^*)^{-1} (r-Aw). \]

%%
% Tolerance and number of iterations for the conjugate gradient.

opts.epsilon = 1e-9; 
opts.niter_max = 150;

%%
% Adapt conjugate gradient fucntion to handle variables that are not
% vectors.

flat = @(x)x(:);
resh = @(x)reshape(x, [n n p+2]);
mycg = @(B,y)resh( perform_cg(@(r)flat(B(resh(r))),y(:),opts) );

%%
% The operator \((A A^*)^{-1}\) can be computed using conjugate gradient.

pA = @(r)mycg(@(s)A(AS(s)),r);

%%
% Define the projection operator \(\text{Prox}_{\la G}\).

ProxG = @(w,lambda)w + AS( pA(r0-A(w)) );

%%
% Check that \(\text{Prox}_{\la G}\) implements the projection on the
% constraint \(Aw=y\).

w = randn(n,n,p,3);
err = @(w)mynorm(A(w)-r0)/mynorm(r0);
fprintf('Error before projection: %.2e\n', err(w));
fprintf('Error before projection: %.2e\n', err(ProxG(w)));

%% Douglas-Rachford Solver
% The discrete optimal transport problem can be written as
% \[ \umin{w} J(w) + G(w) \qwhereq G=\iota_{\Cc}.  \]

%%
% The Douglas-Rachford (DR) algorithm is an iterative scheme to minimize
% functionals of the form \(J+G\)
% where \(J\) and \(G\) are convex functions for which one is able to
% comptue the proximal operators.

%%
% A DR iteration reads
% \[ \tilde w_{k+1} = \pa{1-\frac{\mu}{2}} \tilde w_k + 
%   \frac{\mu}{2} \text{rPox}_{\gamma J}( \text{rProx}_{\gamma G}(\tilde w_k)  ) 
%   \qandq w_{k+1} = \text{Prox}_{\gamma G}(\tilde w_{k+1}). \]

%%
% We have use the following shortcuts:
% \[   \text{rProx}_{\gamma G}(w) = 2\text{Prox}_{\gamma F}(w)-w \]

%%
% One can show that for any value of \(\gamma>0\), any \( 0 < \mu < 2 \), 
% and any \(\tilde w_0\), \(w_k \rightarrow w^\star\)
% which is a minimizer of the minimization of \(J+G\).

%%
% To learn more about this algorithm, you can read:

%%
% _Proximal Splitting Methods in Signal Processing_, Patrick L. Combettes
% and Jean-Christophe Pesquet, in: Fixed-Point Algorithms for Inverse
% Problems in Science and Engineering, New York: Springer-Verlag, 2010.


%%
% Set the value of \(\mu\) and \(\gamma\).
% You might consider using your own values to speed up the convergence.

mu = 1;
gamma = 1;

%%
% Define the rProx operators.

rProxJ = @(w,tau)2*ProxJ(w,tau)-w;
rProxG = @(w,tau)2*ProxG(w,tau)-w;

%%
% Number of iterations.

niter = 200;

%%
% Initialization using linear interpolation of the densities.

t = repmat( reshape(linspace(0,1,p), [1 1 p]), [n n 1]);
f = (1-t) .* repmat(f0, [1 1 p]) + t .* repmat(f1, [1 1 p]);
m = zeros(n,n,p,2);
w0 = cat(4, m,f);

%%
% Display the initialization.

sel = round(linspace(1,p,6));
clf;
imageplot( mat2cell(w0(:,:,sel,3), n, n, ones(6,1)) , '', 2,3);

% 
w = ProxG(w0,gamma);
mynorm(A(w0)-r0)/mynorm(r0)
mynorm(A(w)-r0)/mynorm(r0)
%

%EXO
%% Implement the DR iterative algorithm on |niter| iterations.
%% Keep track of the evolution of the energy \(J\).
energy = [];
constr = [];
tw = w0;
for i=1:niter
    tw_old = tw;
    w = ProxG(tw,gamma);
    rw = 2*w-tw_old;
    tw = (1-mu/2)*tw + mu/2*rProxJ( rw, gamma );
    % 
    energy(i) = J(w);
    constr(i) = mynorm( A(w)-r0 ) / mynorm(r0); 
end
clf;
h = plot(min(energy, energy(1)));
set(h, 'LineWidth', 2);
title('J(w)');
axis tight;
%EXO

%%
% Display the resulting density \(f(x,t)\) for \(t\) from 0 to 1.

sel = round(linspace(1,p,6));
clf;
imageplot( mat2cell(w(:,:,sel,3), n, n, ones(6,1)) , '', 2,3);

