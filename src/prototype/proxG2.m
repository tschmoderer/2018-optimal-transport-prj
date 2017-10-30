function [mbart,fbart,mt,ft] = proxG2(mbar,fbar,m,f)
globals;

mbar1 = m*Interpm_adj;
fbar1 = Interpf_adj*f;


mbart = (mbar + mbar1)\(eye(N+1) + Interpm*Interpm_adj);
fbart = (fbar + fbar1)\(eye(Q+1) + Interpf_adj*Interpf);

mt = mbart*Interpm;
ft = Interpf*fbart;
end