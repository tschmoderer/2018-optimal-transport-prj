function [mbart,fbart,mt,ft] = proxG2(mbar,fbar,m,f)
    globals;

    Am = (eye(N+2) + Interpm*Interpm_adj)';
    Af = eye(Q+2) + Interpf_adj*Interpf;

    Bm = (mbar+m*Interpm_adj)';
    Bf = fbar + Interpf_adj*f;

    fbart = Af\Bf;
    mbart = (Am\Bm)';

    mt = mbart*Interpm;
    ft = Interpf*fbart;
end