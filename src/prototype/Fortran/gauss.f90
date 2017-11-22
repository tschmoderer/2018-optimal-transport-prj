function gauss(mu,sigma,N)
    implicit none
    integer, intent(in) :: N;
    double precision, intent(in) :: mu, sigma;
    double precision, dimension(N+1):: gauss; 
    double precision, dimension(N+1) :: x; 
    integer :: i;

    do i = 0,N 
        x = 1.0*i/(1.0*N);
    end do
    x = (x-mu)/(1.0*sigma);
    x = -0.5*x*x;

    gauss = exp(x);
    return
end function 