% Validée
function [mbart,fbart,mt,ft] = proxG2(mbar,fbar,m,f)
    globals;
    %% Matrice d'interpolation %%
    Interpm = [diag(ones(1,N+1));zeros(1,N+1)] + [zeros(1,N+1);diag(ones(1,N+1))];
    Interpm = Interpm/2; % m = mbar*Im
    Interpf = [diag(ones(1,Q+1)) zeros(Q+1,1)] + [zeros(Q+1,1) diag(ones(1,Q+1))];
    Interpf = Interpf/2;

    Interpm_adj = Interpm';
    Interpf_adj = Interpf';
    
    Am = eye(N+2) + Interpm*Interpm_adj; % sym def positive
    Af = eye(Q+2) + Interpf_adj*Interpf;
    
    Bm = mbar + m*Interpm_adj;
    Bf = fbar + Interpf_adj*f;

    fbart = Af\Bf;
    mbart = Bm/Am;

    mt = mbart*Interpm;
    ft = Interpf*fbart;
end

