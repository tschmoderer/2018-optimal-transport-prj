function [pC, error] = projC(w)
    globals;

    dx  = @(m1)  N*(m1(:,[2:end end],:) - m1); % dérivation selon ---> x
    dy  = @(m2)  P*(m2([2:end end],:,:) - m2); % derivation selon ---> y
    dt  = @(f)   Q*(f(:,:,[2:end end]) - f); % dérivation selon ---> t 

    dxS = @(dm1) N*cat(2, -dm1(:,1,:) , dm1(:,1:end-2,:) - dm1(:,2:end-1,:) , dm1(:,end-1,:));
    dyS = @(dm2) P*cat(1, -dm2(1,:,:) , dm2(1:end-2,:,:) - dm2(2:end-1,:,:) , dm2(end-1,:,:));
    dtS = @(df)  Q*cat(3, -df(:,:,1) , df(:,:,1:end-2) - df(:,:,2:end-1) , df(:,:,end-1));

    %% check adjoint -- ok
%     dxm = dx(m(:,:,:,1)); dym = dy(m(:,:,:,2));
% 
%     rdxm = rand(size(dxm)); rdym = rand(size(dym));
%     dxSrdxm = dxS(rdxm); dySrdym = dyS(rdym);
% 
%     dtf = dt(f); 
%     rdtf = rand(size(dtf));
%     dtSrdtf = dtS(rdtf);
% 
%     sum(sum(sum(m(:,:,:,1).*dxSrdxm - dxm.*rdxm)))
%     sum(sum(sum(m(:,:,:,2).*dySrdym - dym.*rdym)))
%     sum(sum(sum(f.*dtSrdtf - dtf.*rdtf)))

    %% construction opérateurs 
    grad = @(u) cat(4,dx(u),dy(u));
    div  = @(w) -dxS(w(:,:,:,1)) - dyS(w(:,:,:,2));

    A    = @(w)  cat(3,div(w(:,:,:,1:2)) + dt(w(:,:,:,3)), w(:,:,end,3) , w(:,:,1,3));
    U    = @(y0,y1) cat(3,y1,zeros(P+1,N+1,Q-2),y0);
    AS   = @(Aw) cat(4,-grad(Aw(:,:,1:Q)),dtS(Aw(:,:,1:Q)) + U(Aw(:,:,end-1),Aw(:,:,end)));

    %% check adjoint -- ok
%     divW = div(m);
%     rdivW = rand(size(divW));
%     g = grad(rdivW);
% 
%     sum(sum(sum(sum(m.*g)))) + sum(sum(sum(divW.*rdivW))) % un + car grad* = -div /!\
% 
%     Aw = A(w);
% 
%     rAw = rand(size(Aw));
%     ASrAw = AS(rAw);
% 
%     sum(sum(sum(sum(w.*ASrAw)))) - sum(sum(sum(sum(Aw.*rAw))))

    %% opérateurs de projetction
    y = cat(3, zeros(P+1,N+1,Q),f0,f1); % second membre

    flat = @(x) x(:);
    resh = @(x) reshape(x,P+1,N+1,Q+2);

    do_cg =@(B,y) resh(cg(@(r)flat(B(resh(r))),y(:))); % solve B*x = y with CG
    pA = @(r) do_cg(@(s)A(AS(s)),r); % solve (A*A')*x = r

    pC = w + AS(pA(y - A(w)));

    mynorm = @(x) norm(x(:));
    err = @(w) mynorm(A(w)-y)/norm(y(:));
    error = err(pC);

    %% check div=0
%     w = rand(P+1,N+1,Q,3);
% 
%     fprintf('Error before projection: %.2e\n', err(w));
%     fprintf('Error before projection: %.2e\n', err(w + AS(pA(y - A(w)))));

end

