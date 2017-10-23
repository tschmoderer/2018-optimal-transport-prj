function [mbart,fbart,mt,ft] = proxG2(mbar,fbar,m,f)

mbar1 = m*Interpm_adj;
fbar1 = Interpf_adj*f;


mbart = (mbar + mbar1)\(eye(N+1) + Interpm_adj*Interpm);
fbart = (fbar + fbar1)\(eye(N+1) + Interpf_adj*Interpf);

mt = mbart*Interpm;
ft = Interpf*fbart;
end