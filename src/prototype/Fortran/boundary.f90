!! Etant donné mbar et fbar !!
!! Renvoi les fontieres !! 

subroutine boundary(mbar,fbar,mleft,mright,fup,fdown,N,Q)
    implicit none
    integer, intent(in) :: N,Q;
    real, dimension(N+2,Q+1), intent(in) :: mbar;
    real, dimension(N+1,Q+2), intent(in) :: fbar;
    real, dimension(Q+1), intent(out) :: mleft,mright;
    real, dimension(N+1), intent(out) :: fup,fdown;

    mleft = mbar(1,:);
    mright = mbar(N+2,:);
    fup = fbar(:,1);
    fdown = fbar(:,Q+2);
end subroutine