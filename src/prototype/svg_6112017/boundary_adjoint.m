% ordre des fontieres
%       z
% x |       |y
%       t

function [mbar, fbar] = boundary_adjoint(x,y,z,t,N,Q)

mbar = [x , zeros(Q+1,N) , y];
fbar = [z ; zeros(Q,N+1) ; t];

end