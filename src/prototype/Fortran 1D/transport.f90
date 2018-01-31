program transport
use procedures
    implicit none

    real, dimension(1:N+1) :: f0,f1
    real, dimension(1:Q+1,1:N+1,2) :: z = 0, w0 = 0, w1 = 0
    real, dimension(1:niter) :: cout = 0, minF = 0, div = 0

    integer :: i,j

	f0 = normalise(eps + gauss(0.5,0.05))
	f1 = normalise(eps + gauss(0.5,0.05))

    do i = 1,Q+1
        w0(i,:,2) = normalise(w0(i,:,2) + eps) 
        z(i,:,2)  = normalise(z(i,:,2) + eps) 
    end do 

    do i = 1,niter
		w1 = w0 + alpha*(proxJ(2.0*z - w0) - z)
		call projC(z,div(i),w1,f0,f1)

        cout(i) = cost(z)
        minF(i) = minval(z(:,:,2))

		print *, i, sum(w0(:,:,2)), sum(w1(:,:,2)), sum(z(:,:,2)), cout(i), minF(i), div(i)
        
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
            write(1,*) (j-1)/(1.0*N), (i-1)/(1.0*Q), w0(i,j,2)
        end do
    end do
    close(1)

    open(1,file='results/w1.dat');
    write(1,*) "# ", "X ", "T ", "Z "
    do i = 1,Q+1
        do j = 1,N+1
            write(1,*) (j-1)/(1.0*N), (i-1)/(1.0*Q), w1(i,j,2)
        end do
    end do
    close(1)

    open(1,file='results/transport.dat');
    write(1,*) "# ", "X ", "T ", "Z "
    do i = 1,Q+1
        do j = 1,N+1
            write(1,*) (j-1)/(1.0*N), (i-1)/(1.0*Q),  z(i,j,2)
        end do
    end do
    close(1)

end program
