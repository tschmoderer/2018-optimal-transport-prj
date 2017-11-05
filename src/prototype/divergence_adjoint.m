function [mbar, fbar] = divergence_adjoint(d)
globals;
mbar = [-d(:,1) d(:,1:end-1)-d(:,2:end) d(:,end)];
fbar = [-d(1,:) ; d(1:end-1,:)-d(2:end,:) ; d(end,:)];

end