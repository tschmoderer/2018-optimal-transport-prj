!function [mbar fbar] = interpolation_adjoint(m,f)
!mbar = [m(:,1) m(:,2:end)+m(:,1:end-1) m(:,end)];
!fbar = [f(1,:) ; f(2:end,:)+f(1:end-1,:) ; f(end,:)];

!mbar = 0.5*mbar;
!fbar = 0.5*fbar;
!end

subroutine interpolation_adjoint(m,f,mbar,fbar,N,Q)
    integer, intent(in) :: N,Q;
    real, dimension(Q+1,N+2), intent(out) :: mbar;
    real, dimension(Q+2,N+1), intent(out) :: fbar;
    real, dimension(Q+1,N+1), intent(in) :: m,f;

    mbar(:,1) = m(:,1); mbar(:,N+2) = m(:,N+1);
    mbar(:,2:N+1) = m(:,2:N+1) + m(:,1:N); 
   
    fbar(1,:) = f(1,:); fbar(Q+2,:) = f(Q+1,:);
    fbar(2:Q+1,:) = f(2:Q+1,:) + f(1:Q,:);

    mbar = 0.5*m;
    fbar = 0.5*f;
end subroutine