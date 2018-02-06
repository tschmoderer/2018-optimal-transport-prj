% Fonction projC
% Out : 
%     - pC : la projection sur les contraintes de w
%      - error : l'erreur commise sur le respect des contraintes
% In  : 
%     - w  : le point d'évaluation 
% Timothée Schmoderer 
% INSA Rouen Normandie 2017/2018

function [pC, error] = projC(w)
    globals;
    dx  = @(m) N*(m(:,[2:end end]) - m); % dérivation selon ---> x
    dt  = @(f) Q*(f([2:end end],:) - f); % dérivation selon ---> t 
    dxS = @(dm) N*[-dm(:,1) , dm(:,1:end-2) - dm(:,2:end-1) , dm(:,end-1)];
    dtS = @(df) Q*[-df(1,:) ; df(1:end-2,:) - df(2:end-1,:) ; df(end-1,:)];
    
    div  = @(w) -dxS(w(:,:,1)) + dt(w(:,:,2));
    grad = @(d) cat(3,-dx(d(:,:)),dtS(d(:,:)));

    A  = @(w)  [-dxS(w(:,:,1)) + dt(w(:,:,2)); w(end,:,2) ;w(1,:,2)];
    AS = @(Aw) cat(3,-dx(Aw(1:Q+1,:)),dtS(Aw(1:Q+1,:)) + [Aw(Q+3,:);zeros(Q-1,N+1);Aw(Q+2,:)]);

    y = [zeros(Q+1,N+1);f0;f1];

%     % gradient conjugué
%     b = y - A(w);
%     x = zeros(size(b));
%     r = b - A(AS(x));
%     p = r;
%     rsold = sum(r(:).*r(:));
%     for i = 1:prod(size(b))
%         Ap = A(AS(p));
%         alpha = rsold / sum(p(:).*Ap(:));
%         x = x + alpha * p;
%         r = r - alpha * Ap;
%         rsnew = sum(r(:).*r(:));
%         if sqrt(rsnew) < 1e-10
%               break;
%         end
%         p = r + (rsnew / rsold) * p;
%         rsold = rsnew;
%     end
    
    b = div(w);
    u = poisson2d_Neumann(-b(2:end-1,2:end-1)); 
    gu = grad(u); gu0 = gu; gu0(2:end-1,2:end-1) = 0;
    pC = zeros(size(w));
   
    pC(2:end-1,2:end-1,:) = w(2:end-1,2:end-1,:) + gu - gu0;
    
    pC(2:end-1,:,1) = w(2:end-1,:,1) - diff(u,[],1);
    pC(:,2:end-1,2) = w(:,2:end-1,2) - diff(u,[],2);
    
    pC(:,[1 end],1) = 0; pC([1 end],:,2) = [f1 ; f0];

    error = norm(A(pC) - y)/norm(y);
end

