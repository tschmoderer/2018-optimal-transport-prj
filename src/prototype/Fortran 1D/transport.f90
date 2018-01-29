program transport
use procedures
    implicit none
    integer, parameter :: N = 21, Q = 19, niter = 94
    double precision, parameter :: epsilon = 1e-10, alpha = 1.0, gamma = 1.0
    double precision, dimension(1:N+1) :: f0,f1
    double precision, dimension(1:Q+1,1:N+1,2) :: z = 0, w0 = 0, w1 = 0
    double precision, dimension(1:niter) :: cout = 0, minF = 0, div = 0

    integer :: i,j,k

	f0 = normalise(epsilon + gauss(0.5d0,0.05d0,N),N)
	f1 = normalise(epsilon + gauss(0.5d0,0.05d0,N),N)

    do i = 1,niter
		w1 = w0 + alpha*(proxJ(2.0*z-w0,gamma,N,Q)-z)
		call projC(z,div(i),w1,f0,f1,N,Q)

        cout(i) = cost(z,N,Q)
        minF(i) = minval(z(:,:,2))
		if (modulo(i,1).EQ. 0) then
			print *, i, ' ',sum(w0), sum(w1), sum(z), cout(i), minF(i), div(i)
		end if
        
        w0 = w1;
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