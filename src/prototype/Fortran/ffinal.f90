subroutine ffinal(x,N,f1)
    implicit none
    integer, intent(in) :: N;
    double precision, dimension(N+1), intent(in) :: x;
    double precision, dimension(N+1), intent(out) :: f1;
    double precision, parameter :: mu = 0.9, s = 0.005, minimal = 0.000001;
    
    f1 = minimal + exp(-0.5*((x-mu)/(1.0*s))**2);
    f1 = f1/(1.0*sum(f1));
end subroutine