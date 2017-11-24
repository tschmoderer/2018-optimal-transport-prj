!function gauss(mu,sigma,N)
subroutine gauss(mu,sigma,N,G)
    implicit none
    integer, intent(in) :: N;
    double precision, intent(in) :: mu, sigma;
    double precision, dimension(N+1), intent(out) :: G; 
    double precision, dimension(N+1) :: x; 
    integer :: i;

    do i = 0,N 
        x = 1.0*i/(1.0*N);
    end do
    x = (x-mu)/(1.0*sigma);
    x = -0.5*x*x;

    G = dexp(x);
end subroutine
!end function 
