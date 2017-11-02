subroutine finitial(x,N,f0)
    implicit none
    integer, intent(in) :: N;
    double precision, dimension(N+1), intent(in) :: x;
    double precision, dimension(N+1), intent(out) :: f0;
    double precision, parameter :: mu = 0.1, s = 0.005, minimal = 0.000001;
    
    f0 = minimal + exp(-0.5*((x-mu)/(1.0*s))**2);
    f0 = f0/(1.0*sum(f0)); ! Normalise
end subroutine