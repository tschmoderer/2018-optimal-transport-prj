subroutine finitial(x,N,f0)
    implicit none
    integer, intent(in) :: N;
    real, dimension(N+1), intent(in) :: x;
    real, dimension(N+1), intent(out) :: f0;
    real :: mu = 0.1;
    real :: s = 0.005;
    real :: pi=4.D0*DATAN(1.D0)
    
    f0 = exp(-0.5*((x-mu)/(1.0*s))**2)/(1.0*s*sqrt(2.0*pi));
end subroutine