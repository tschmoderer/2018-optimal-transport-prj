subroutine normalise(f,res,N)
    integer, intent(in) :: N
    double precision, dimension(1:N+1), intent(in) :: f
    double precision, dimension(1:N+1), intent(out) :: res
    res = f/sum(f)
end subroutine

subroutine gauss(f,mu,sigma,N)
    real, intent(in) :: mu, sigma
    integer, intent(in) :: N 
    double precision, dimension(1:N+1), intent(out) :: f
    integer :: i
    do i = 1,N+1 
        f(i) = exp(-0.5*(((i-1)/(1.0*N) - mu)/sigma)**2)
    end do
end subroutine
    
subroutine cost(c,w,N,Q)
    integer, intent(in) :: N,Q
    double precision, dimension(1:Q+1,1:N+1,2), intent(in) :: w
    double precision, intent(out) :: c
    c = sum(w(:,:,1)**2/w(:,:,2))
end subroutine

subroutine proxJ(Pw,w,g,N,Q)
    integer, intent(in) :: N,Q
    double precision, intent(in) :: g
    double precision, dimension(1:Q+1,1:N+1,2), intent(in) :: w
    double precision, dimension(1:Q+1,1:N+1,2), intent(out) :: Pw

    double precision, dimension(1:Q+1,1:N+1) :: mt, ft
    double precision, dimension(1:Q+1,1:N+1) :: x0, x1, poly, dpoly

    integer :: k = 0, i, j
    x0 = 1000; x1 = 2000
    mt = w(:,:,1); ft = w(:,:,2)
    ! Newton
    do while (maxval(dabs(x0-x1)) > 1e-5 .AND. k < 1500)
        x0 = x1;
        poly  = (x0 - ft)*((x0 + g)**2)-0.5*g*(mt**2)
        dpoly = 2*(x0 + g)*(x0 - ft) + (x0 + g)**2
        x1 = x0 - poly/dpoly
        k = k+1;
    end do 
    
    do i = 1,Q+1
        do j = 1,N+1 
            if (x1(i,j) .GT. 0) then
                Pw(i,j,2) = x1(i,j)
                Pw(i,j,1) = Pw(i,j,2)*mt(i,j)/(Pw(i,j,2) + g)
            else 
                Pw(i,j,2) = 0
                Pw(i,j,1) = 0
            end if
        end do 
    end do
end subroutine

program transport
    implicit none
    integer, parameter :: N = 21, Q = 19, niter = 50
    double precision, parameter :: epsilon = 1e-10, alpha = 1.0, beta = 1.0, gamma = 1.0
    double precision, dimension(1:N+1) :: f0,f1
    double precision, dimension(1:Q+1,1:N+1,2) :: z = 0, w0 = 0, w1 = 0
    double precision, dimension(1:niter) :: cout = 0, minF = 0, divV = 0

    integer :: i,j

    call gauss(f0,0.5,0.05,N)
    call gauss(f1,0.8,0.05,N)

    call normalise(epsilon + f0,f0,N)
    call normalise(epsilon + f1,f1,N)

    do i = 1,niter
        call proxJ(w1,2*z - w0,gamma,N,Q)
        w1 = w0 + alpha*(w1 - z)
        call projC(z,divV(i),w1,f0,f1,N,Q)
        print *, 'iteration :', i, ' ',sum(w0), sum(w1), sum(z)
        
        w0 = w1;


        call cost(cout(i),z,N,Q)
        minF(i) = minval(z(:,:,2))
    end do

    open(1,file='results/f0.dat');
    write(1,*) "# ", "X ", "Y "
    do i = 1,N+1
        write(1,*) (i-1)/(1.0*N), f0(i)
    end do
    close(1)

    open(1,file='results/f1.dat');
    write(1,*) "# ", "X ", "Y "
    do i = 1,N+1
        write(1,*) (i-1)/(1.0*N), f1(i)
    end do
    close(1)


    open(1,file='results/w0.dat');
    write(1,*) "# ", "X ", "T ", "Z "
    do i = 1,Q+1
        do j = 1,N+1
            write(1,*) (i-1)/(1.0*Q), (j-1)/(1.0*N), w0(i,j,2)
        end do
    end do
    close(1)

    open(1,file='results/w1.dat');
    write(1,*) "# ", "X ", "T ", "Z "
    do i = 1,Q+1
        do j = 1,N+1
            write(1,*) (i-1)/(1.0*Q), (j-1)/(1.0*N), w1(i,j,2)
        end do
    end do
    close(1)

    open(1,file='results/transport.dat');
    write(1,*) "# ", "X ", "T ", "Z "
    do i = 1,Q+1
        do j = 1,N+1
            write(1,*) (i-1)/(1.0*Q), (j-1)/(1.0*N), z(i,j,2)
        end do
    end do
    close(1)

end program