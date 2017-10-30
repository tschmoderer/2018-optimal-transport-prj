!! Etant donné mbar et fbar !!
!! Renvoi les fontieres !! 

subroutine boundary(mbar,fbar,mleft,mright,fup,fdown,N,Q)
    implicit none
    integer, intent(in) :: N,Q;
    double precision, dimension(N+2,Q+1), intent(in) :: mbar;
    double precision, dimension(N+1,Q+2), intent(in) :: fbar;
    double precision, dimension(Q+1), intent(out) :: mleft,mright;
    double precision, dimension(N+1), intent(out) :: fup,fdown;

    mleft = mbar(1,:);
    mright = mbar(N+2,:);
    fup = fbar(:,1);
    fdown = fbar(:,Q+2);
end subroutine