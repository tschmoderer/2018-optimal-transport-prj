program transport
    implicit none
    integer, parameter :: N = 51, Q = 49, niter = 2000
    double precision, parameter :: eps = 1e-10, alpha = 1.0, g = 1.0, b = 1
    double precision, dimension(N+1) :: f0, f1
    double precision, dimension(Q+1,N+1,2) :: z = 0, w0 = 0, w1 = 0
	integer, dimension(Q+1,N+1) :: obstacle = 0
    double precision, dimension(niter) :: cout, minF
    integer :: i,k
	 
    f0 = normalise(eps + gauss(0.2d0,0.05d0))
    f1 = normalise(eps + gauss(0.8d0,0.05d0))

	obstacle(25,:) = 1; obstacle(25,2:4) = 0
!	obstacle(7,:)  = 1; obstacle(7,12:14) = 0
!	obstacle(10,:) = 1; obstacle(10,2:4) = 0
!	obstacle(4,:) = 1; obstacle(4,17:19) = 0

    do i = 1,niter
        w1 = w0 + alpha*(proxJ(2*z-w0) - z)
        z = projC(w1) 
		w0 = w1
		
        cout(i) = J(z)
       
        if (modulo(i,50) .EQ. 0) print *, i, cout(i)
        minF(i) = minval(z(:,:,2))
    end do 

    open(1,file='results/transport.dat');
    write(1,*) "# ", "X ", "T ", "Z "
    do i = 1,Q+1
        do k = 1,N+1
           write(1,*) (k-1)/(1.0*N), (Q - i + 1)/(1.0*Q),  z(i,k,2)
        end do
    end do
    close(1)
    
    open(1,file='results/energie.dat');
    write(1,*) "# ", "X ", "T ", "Z "
    do i = 1,niter
		write(1,*) i, cout(i)
    end do
    close(1)    
    
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

	open(8,file="results/plot.gnu"); 
	write(8,*) 'set contour' 
	write(8,*) 'set cntrparam levels 30'
	write(8,*) 'unset key'
	write(8,*) 'set pm3d'
	write(8,*) 'unset colorbox'
	write(8,*) 'set hidden3d'
	write(8,*) 'set title "Transport Optimal"'
	write(8,*) 'set xlabel "x"'
	write(8,*) 'set ylabel "t"'
	write(8,*) 'set dgrid3d ', Q+1, ',', N+1
	write(8,*) 'splot "results/transport.dat" with lines'
	close(8);
	call system("gnuplot -p results/plot.gnu");
    
    contains

!! Gauss
    function gauss(mu,sigma) result(f) 
	implicit none
        double precision :: mu, sigma 
        double precision, dimension(N+1) :: f
        integer :: i
        do i = 1,N+1
            f(i) = exp(-0.5*((((i-1)/(1.0*N))-mu)/sigma)**2)
        end do 
    end function gauss

!! Normalise
    function normalise(f) result(nf) 
	implicit none
	    double precision, dimension(N+1) :: f, nf
        nf = f/sum(f)
    end function

!! Cout 
    function J(w) result(c) 
    implicit none
        double precision, dimension(Q+1,N+1,2) :: w
        double precision :: c
        c = 0.5*sum(w(:,:,1)**2/max(w(:,:,2),eps,1e-10)**b)
    end function J 

!! Proximal de J 
    function proxJ(w) result(pw)
    implicit none
        double precision, dimension(Q+1,N+1,2) :: w, pw
        double precision, dimension(Q+1,N+1) :: mt, ft, x0, x1, poly, dpoly
        integer :: k
        
        mt = w(:,:,1); ft = w(:,:,2);
        x0 = 1; x1 = 2; k = 0;

        do while (maxval(dabs(x0-x1)) .GT. 1e-5  .AND. k .LT. 1500)
            x0 = x1
            if (b .EQ. 1) then ! Cas transport
				poly  = (x0-ft)*(x0+g)**2 - 0.5*g*mt**2
				dpoly = 2*(x0+g)*(x0-ft) + (x0+g)**2
			else if (b .EQ. 0) then ! Interpolation L2
				x1 = ft
				exit
			else 
				poly = x0**(1.0-b)*(x0-ft)*((x0**b+g)**2)-0.5*b*g*mt**2
				dpoly = (1.0-b)*x0**(-b)*(x0-ft)*((x0**b+g)**2) + x0**(1-b)*((x0**b+g)**2 +2*b*(x0-ft)*x0**(b-1)*(x0**b+g) )
			end if
			
			where (x0 .GT. eps) x1 = x0 - poly/dpoly
			where (x0 .LT. eps) x1 = eps		
			
            k = k+1
        end do

        where (x1 .LT. eps) x1 = eps
        where (obstacle .GT. 0) x1 = eps
        pw(:,:,2) = x1
        pw(:,:,1) = (x1**b)*mt/(x1**b+g) 
    end function proxJ

!! Projection sur C 
    function projC(w) result(pc)
    implicit none
        double precision, dimension(Q+1,N+1,2) :: w, pc
        double precision, dimension(Q+3,N+1) :: y, b, x, r, p, Ap
        double precision :: alpha, rold, rnew
        integer :: i 
        
        y = 0
        y(Q+2,:) = f0
        y(Q+3,:) = f1
        
     !! Gradient conjugué
        x = 0
        b = y - A(w)
        r = b - A(AS(x))
        p = r
        rold = sum(r*r)
        do i = 1,(N+1)*(Q+3)
			Ap = A(AS(p))
			alpha = rold/sum(p*Ap)
			x = x + alpha*p
			r = r - alpha*Ap
			rnew = sum(r*r)
			if (dsqrt(rnew) .LT. 1e-10) exit
			p = r + (rnew/rold)*p
			rold = rnew
        end do
        pc = w + AS(x) 
    end function projC
    
!! Dérivation en x
	function dx(m) result(dm)
	implicit none
		double precision, dimension(Q+1,N+1) :: m, dm
		dm = 0
		dm(:,1:N) = N*(m(:,2:N+1) - m(:,1:N))
	end function dx
	
!! Dérivation en t
	function dt(f) result(df)
	implicit none
		double precision, dimension(Q+1,N+1) :: f, df
		df = 0
		df(1:Q,:) = Q*(f(2:Q+1,:) - f(1:Q,:))
	end function dt

!! Adjoint de la dérivation en x 
	function dxS(dm) result(m)
	implicit none
		double precision, dimension(Q+1,N+1) :: m, dm
		m(:,1) = -N*dm(:,1)
		m(:,2:N) = N*(dm(:,1:N-1) - dm(:,2:N))
		m(:,N+1) = N*dm(:,N)		
	end function dxS

!! Adjoint de la dérivation en t 
	function dtS(df) result(f)
	implicit none
		double precision, dimension(Q+1,N+1) :: f, df
		f(1,:) = -Q*df(1,:)
		f(2:Q,:) = Q*(df(1:Q-1,:) - df(2:Q,:))
		f(Q+1,:) = Q*df(Q,:)
	end function dtS		

!! Opérateur A
	function A(w) result(aw)
	implicit none
		double precision, dimension(Q+1,N+1,2) :: w
		double precision, dimension(Q+3,N+1) :: aw
		aw(1:Q+1,:) = -dxS(w(:,:,1)) + dt(w(:,:,2))
		aw(Q+2,:) = w(Q+1,:,2)
		aw(Q+3,:) = w(1,:,2)
	end function A

!! Adjoint de A
	function AS(aw) result(w)
	implicit none
		double precision, dimension(Q+1,N+1,2) :: w
		double precision, dimension(Q+3,N+1) :: aw		
		w(:,:,1) = -dx(aw(1:Q+1,:))
		w(:,:,2) = dtS(aw(1:Q+1,:))
		w(1,:,2) = w(1,:,2) + aw(Q+3,:)
		w(Q+1,:,2) = w(Q+1,:,2) + aw(Q+2,:)
	end function AS
end program transport
