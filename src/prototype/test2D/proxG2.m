function [Ut,Vt] = proxG2(U,V)
    globals;
    
    Ut = pG2*(U + Interp'*V);
    Vt = Interp*Ut;
end

