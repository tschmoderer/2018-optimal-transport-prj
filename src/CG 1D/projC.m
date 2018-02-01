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

    A  = @(w)  [-dxS(w(:,:,1)) + dt(w(:,:,2)); w(end,:,2) ;w(1,:,2)];
    AS = @(Aw) cat(3,-dx(Aw(1:Q+1,:)),dtS(Aw(1:Q+1,:)) + [Aw(Q+3,:);zeros(Q-1,N+1);Aw(Q+2,:)]);

    y = [zeros(Q+1,N+1);f0;f1];

    flat = @(x) reshape(x,(Q+3)*(N+1),1);
    resh = @(x) reshape(x,Q+3,N+1);

    pC = w + AS(resh(cg(@(s)flat(A(AS(resh(s)))),flat(y-A(w)))));

    err = @(w) norm(A(w)-y)/norm(y);
    error = err(pC);
end

