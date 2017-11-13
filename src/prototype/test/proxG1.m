function [pU,pV] = proxG1(U,V,gamma)
    pV = proxJ(V,gamma);
    pU = projC(U);
end

