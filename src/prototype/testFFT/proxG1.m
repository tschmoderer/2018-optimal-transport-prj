function [pU,pV] = proxG1(U,V)
    pV = proxJ(V);
    pU = projC(U);
end

