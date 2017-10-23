function [m f] = boundary_adjoint(x,y,z,t,N,Q)

m = [x zeros(Q+1,N) y];
f = [z ; zeros(Q,N+1) ; t];

end