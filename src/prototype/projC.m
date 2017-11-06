function [Umbar, Ufbar] = projC(mbar,fbar)
globals;

d = divergence(mbar,fbar);
b = boundary(mbar,fbar); 

solu = poisson(d,b.m(:,1),b.m(:,2),b.f(1,:),b.f(2,:),N,Q,1e-3);

[solmbar, solfbar] = divergence_adjoint(solu);

Umbar = mbar - solmbar + Cstmbar;
Ufbar = fbar - solfbar + Cstfbar;

Ufbar(1,:) = fbar(1,:) - solu(1,:) + Cst(1,:); 
Ufbar(end,:) = fbar(end,:) - solu(end,:) + Cst(end,:); 
Umbar(:,1) = mbar(:,1) - solu(:,1) + Cst(:,1);
Umbar(:,end) = mbar(:,end) - solu(:,end) + Cst(:,end);

%  [X Y] = meshgrid(linspace(0,1,N+1),linspace(0,1,Q+2));
%  surf(X,Y,Ufbar);
%  title('Ufbar du pbm de poisson');
%  pause
end