subroutine interpolation(mbar,fbar,m,f,N,Q)

    integer, intent(in) :: N,Q;
    real, dimension(Q+1,N+2), intent(in) :: mbar;
    real, dimension(Q+2,N+1), intent(in) :: fbar;
    real, dimension(Q+1,N+1), intent(out) :: m,f;

    m = mbar(:,2:N+2) + mbar(:,1:N+1);
    f = fbar(2:Q+2,:) + fbar(1:Q+1,:);

    m = 0.5*m;
    f = 0.5*f;
end subroutine