function [Ut,Vt] = proxG2(U,V,N,Q)
    mbar = reshape(U(1:(Q+1)*(N+2)),Q+1,N+2);
    fbar = reshape(U((Q+1)*(N+2)+1:end),Q+2,N+1);
    m = reshape(V(1:(N+1)*(Q+1)),Q+1,N+1);
    f = reshape(V((N+1)*(Q+1)+1:end),Q+1,N+1);
    
    %% Matrice d'interpolation %%
    Interpm = [diag(ones(1,N+1));zeros(1,N+1)] + [zeros(1,N+1);diag(ones(1,N+1))];
    Interpm = Interpm/2; % m = mbar*Im
    Interpf = [diag(ones(1,Q+1)) zeros(Q+1,1)] + [zeros(Q+1,1) diag(ones(1,Q+1))];
    Interpf = Interpf/2;

    Interpm_adj = Interpm';
    Interpf_adj = Interpf';
    Am = (eye(N+2) + Interpm*Interpm_adj)';
    Af = eye(Q+2) + Interpf_adj*Interpf;

    Bm = (mbar+m*Interpm_adj)';
    Bf = fbar + Interpf_adj*f;

    fbart = Af\Bf;
    mbart = (Am\Bm)';

    mt = mbart*Interpm;
    ft = Interpf*fbart;
    
    Vt = [reshape(mt,(N+1)*(Q+1),1);reshape(ft,(N+1)*(Q+1),1)];
    Ut = [reshape(mbart,(N+2)*(Q+1),1);reshape(fbart,(N+1)*(Q+2),1)];
end

