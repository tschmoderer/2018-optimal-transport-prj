function [mbar fbar] = interpolation_adjoint(m,f)
%mbar = [m(:,1) m(:,2:end)+m(:,1:end-1) m(:,end)];
%fbar = [f(1,:) ; f(2:end,:)+f(1:end-1,:) ; f(end,:)];
%
%mbar = 0.5*mbar;
%fbar = 0.5*fbar;
globals; 
mbar = m*Interpm_adj;
fbar = Interpf_adj*f;
end