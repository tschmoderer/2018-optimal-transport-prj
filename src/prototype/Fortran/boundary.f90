!! Etant donné mbar et fbar !!
!! Renvoi les fontieres !!
!! Timothée Schmoderer !! 
!! cc 2017!!

subroutine boundary(mbar,fbar,mleft,mright,fup,fdown,N,Q)
    implicit none
    integer, intent(in) :: N,Q;
    double precision, dimension(Q+1,N+2), intent(in) :: mbar;
    double precision, dimension(Q+2,N+1), intent(in) :: fbar;
    double precision, dimension(Q+1), intent(out) :: mleft,mright;
    double precision, dimension(N+1), intent(out) :: fup,fdown;

    mleft = mbar(:,1); ! la première colonne
    mright = mbar(:,N+2); ! la dernière colonne
    fup = fbar(1,:); ! la première ligne 
    fdown = fbar(Q+2,:); ! la dernière ligne
end subroutine