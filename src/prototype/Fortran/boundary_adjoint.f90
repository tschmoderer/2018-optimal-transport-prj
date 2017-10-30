!! Etant donné des fontieres !!
!! reconstruit des variables décentrées !!

subroutine boundary_adjoint(mbar,fbar,mleft,mright,fup,fdown,N,Q)
    implicit none
    integer, intent(in) :: N,Q;
    double precision, dimension(N+2,Q+1), intent(out) :: mbar;
    double precision, dimension(N+1,Q+2), intent(out) :: fbar;
    double precision, dimension(Q+1), intent(in) :: mleft,mright;
    double precision, dimension(N+1), intent(in) :: fup,fdown;

    mbar = 0;
    fbar = 0;
 
    mbar(1,:) = mleft;
    mbar(N+2,:) = mright;

    fbar(:,1) = fup;
    fbar(:,Q+2) = fdown;
end subroutine