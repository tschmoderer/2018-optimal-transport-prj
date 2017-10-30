!! Etant donnée des valeurs décentrées bar et fbar
!! Calcul leur divergence d
subroutine divergence(mbar,fbar,d,N,Q)
    implicit none
    integer, intent(in) :: N,Q;
    double precision, dimension(N+2,Q+1), intent(in) :: mbar;
    double precision, dimension(N+1,Q+2), intent(in) :: fbar;
    double precision, dimension(N+1,Q+1), intent(out) :: d; 

    d = mbar(2:N+2,:) - mbar(1:N+1,:) + fbar(:,2:Q+2) - fbar(:,1:Q+1);
end subroutine