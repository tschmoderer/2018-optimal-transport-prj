function [pmbar,pfbar,pm,pf] = proxG1(mbar,fbar,m,f)
    [pm,pf] = proxJ(m,f);
    [pmbar,pfbar] = projC(mbar,fbar);
end

