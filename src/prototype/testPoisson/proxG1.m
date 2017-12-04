function [pmbar,pfbar,pm,pf] = proxG1(mbar,fbar,m,f,N,Q)
    [pm,pf] = proxJ(m,f);
    [pmbar,pfbar] = projC(mbar,fbar,N,Q);
end

