!! Etant donné des fontieres !!
!! reconstruit des variables décentrées !!

subroutine boundary_adjoint(mbar,fbar,mleft,mright,fup,fdown,N,Q)
    implicit none
    integer, intent(in) :: N,Q;
    double precision, dimension(Q+1,N+2), intent(out) :: mbar;
    double precision, dimension(Q+2,N+1), intent(out) :: fbar;
    double precision, dimension(Q+1), intent(in) :: mleft,mright;
    double precision, dimension(N+1), intent(in) :: fup,fdown;

    mbar = 0;
    fbar = 0;
 
    mbar(:,1) = mleft; ! la première colonne 
    mbar(:,N+2) = mright; ! la dernière colonne 

    fbar(1,:) = fup; !  la première ligne 
    fbar(Q+2,:) = fdown; ! la dernière ligne 
end subroutine