program transport
    implicit none
    integer, parameter :: N = 30, P = 30, Q = 50, niter = 1000
    double precision, parameter :: eps = 1e-10, alpha = 1.98, g = 1./230, b = 0.5
    double precision, dimension(P+1,N+1) :: f0, f1
    double precision, dimension(P+1,N+1,Q,3) :: z = 0, w0 = 0, w1 = 0
    double precision, dimension(niter) :: cout, minF
    integer :: i,k,l
		character(10) :: charI;

	open(1,file='input/f0.dat') 
	do i = 1,P+1
		read(1,*) f0(i,:)
	end do 
	close(1)
	
	open(1,file='input/f1.dat') 
	do i = 1,P+1
		read(1,*) f1(i,:)
	end do 
	close(1)
	
	f0 = normalise(eps + f0)
	f1 = normalise(eps + f1)
	
    f0 = normalise(eps + gauss(0.2d0,0.2d0,0.05d0)) ! + gauss(0.8d0,0.5d0,0.05d0))
    f1 = normalise(eps + gauss(0.8d0,0.8d0,0.05d0)) ! + gauss(0.5d0,0.8d0,0.05d0))

    do i = 1,niter
        w1 = w0 + alpha*(proxJ(2*z-w0) - z)
        z = projC(w1) 
		w0 = w1
		
        cout(i) = J(z)
       
        if (modulo(i,10) .EQ. 0) print *, i, cout(i)
        minF(i) = minval(z(:,:,:,3))
    end do 
    
    open(1,file='results/transport.dat');
    write(1,*) "# ", "X ", "Y ", "T ", "Z "
    do i = 1,P+1 ! y direction
        do k = 1,N+1 ! x direction 
			do l = 1,Q ! t direction
           write(1,*) (k-1)/(1.0*N), (i-1)/(1.0*p), (l-1)/(1.0*Q -1),  z(i,k,l,3)
           end do
        end do
    end do
    close(1)

	do l = 1,Q 
		write(charI,'(I5.5)') l
		open(1,file='results/Transport/'//trim(charI)//'.dat'); 
		write(1,*) "# ", "X ", "Y ", "T ", "Z "
		do i = 1,P+1 ! y direction
			do k = 1,N+1 ! x direction 
				write(1,*) (k-1)/(1.0*N), (i-1)/(1.0*p),  z(i,k,l,3)
			end do
		end do
		close(1)
	end do 

    
    open(1,file='results/energie.dat');
    write(1,*) "# ", "X ", "E "
    do i = 1,niter
		write(1,*) i, cout(i)
    end do
    close(1)    
    
    open(1,file='results/f0.dat');
    write(1,*) "# ", "X ", "Y ", "Z "
    do i = 1,N+1
		do k = 1,P+1
			write(1,*) (i-1)/(1.0*N), (k-1)/(1.0*P),f0(i,k)
        end do 
    end do
    close(1)

    open(1,file='results/f1.dat');
    write(1,*) "# ", "X ", "Y ", "Z "
    do i = 1,N+1
		do k = 1,P+1
			write(1,*) (i-1)/(1.0*N), (k-1)/(1.0*P),f1(i,k)
        end do 
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
	write(8,*) 'set palette gray'
 	write(8,*) 'set view 0,0'
	write(8,*) 'set dgrid3d ', P+1, ',', N+1
	close(8);
!	call system("gnuplot -p results/plot.gnu");    
    
    contains

!! Gauss
    function gauss(muX,muY,sigma) result(f) 
	implicit none
        double precision :: muX, muY, sigma 
        double precision, dimension(P+1,N+1) :: f
        double precision :: tmp
        integer :: i,k
        do i = 1,P+1
			do k = 1,N+1
				tmp = ((k-1)/(1.0*N)-muX)**2 + ((i-1)/(1.0*P)-muY)**2
				f(i,k) = exp(-0.5*tmp/(sigma**2))
			end do
        end do 
    end function gauss

!! Normalise
    function normalise(f) result(nf) 
	implicit none
	    double precision, dimension(P+1,N+1) :: f, nf
        nf = f/sum(f)
    end function

!! Cout 
    function J(w) result(c) 
    implicit none
        double precision, dimension(P+1,N+1,Q,3) :: w
        double precision :: tmp, c
        tmp = sum(w(:,:,:,1)**2 + w(:,:,:,2)**2)
        c = 0.5*sum(tmp/max(w(:,:,:,3),eps,1e-10)**b)
    end function J 

!! Proximal de J 
    function proxJ(w) result(pw)
    implicit none
        double precision, dimension(P+1,N+1,Q,3) :: w, pw
        double precision, dimension(P+1,N+1,Q,2) :: mt
        double precision, dimension(P+1,N+1,Q) :: ft, x0, x1, poly, dpoly
        integer :: k
        
        mt = w(:,:,:,1:2); ft = w(:,:,:,3);
        x0 = 1; x1 = 2; k = 0;

        do while (maxval(dabs(x0-x1)) .GT. 1e-5  .AND. k .LT. 1500)
            x0 = x1
            if (b .EQ. 1) then ! Cas transport
				poly  = (x0-ft)*(x0+g)**2 - 0.5*g*(mt(:,:,:,1)**2 + mt(:,:,:,2)**2)
				dpoly = 2*(x0+g)*(x0-ft) + (x0+g)**2
			else if (b .EQ. 0) then ! Interpolation L2
				x1 = ft
				exit
			else 
				poly = x0**(1.0-b)*(x0-ft)*((x0**b+g)**2)-0.5*b*g*(mt(:,:,:,1)**2 + mt(:,:,:,2)**2)
				dpoly = (1.0-b)*x0**(-b)*(x0-ft)*((x0**b+g)**2) + x0**(1-b)*((x0**b+g)**2 +2*b*(x0-ft)*x0**(b-1)*(x0**b+g) )
			end if
			
			where (x0 .GT. eps) x1 = x0 - poly/dpoly
			where (x0 .LT. eps) x1 = eps		
			
            k = k+1
        end do

        where (x1 .LT. eps) x1 = eps
        
        pw(:,:,:,3) = x1
        pw(:,:,:,1) = (x1**b)*mt(:,:,:,1)/(x1**b+g) 
        pw(:,:,:,2) = (x1**b)*mt(:,:,:,2)/(x1**b+g) 
    end function proxJ

!! Projection sur C 
    function projC(w) result(pc)
    implicit none
        double precision, dimension(P+1,N+1,Q,3) :: w, pc
        double precision, dimension(P+1,N+1,Q+2) :: y, b, x, r, dir, Adir
        double precision :: alpha, rold, rnew
        integer :: i
        y = 0
        y(:,:,Q+1) = f0
        y(:,:,Q+2) = f1
        !! Gradient conjugué
        x = 0
        b = y - A(w)
        r = b - A(AS(x))
        dir = r
        rold = sum(r*r)
        do i = 1,(P+1)*(N+1)*(Q+2)
			Adir = A(AS(dir))
			alpha = rold/sum(dir*Adir)
			x = x + alpha*dir
			r = r - alpha*Adir
			rnew = sum(r*r)
			if (dsqrt(rnew) .LT. 1e-10) exit
			dir = r + (rnew/rold)*dir
			rold = rnew
        end do
        pc = w + AS(x) 
    end function projC
    
!! Dérivation en x
	function dx(m) result(dm)
	implicit none
		double precision, dimension(P+1,N+1,Q) :: m, dm
		dm = 0
		dm(:,1:N,:) = N*(m(:,2:N+1,:) - m(:,1:N,:))
	end function dx
	
!! Dérivation en y
	function dy(m) result(dm)
	implicit none
		double precision, dimension(P+1,N+1,Q) :: m, dm
		dm = 0
		dm(1:P,:,:) = P*(m(2:P+1,:,:) - m(1:P,:,:))
	end function dy
	
!! Dérivation en t
	function dt(f) result(df)
	implicit none
		double precision, dimension(P+1,N+1,Q) :: f, df
		df = 0
		df(:,:,1:Q-1) = Q*(f(:,:,2:Q) - f(:,:,1:Q-1))
	end function dt

!! Adjoint de la dérivation en x 
	function dxS(dm) result(m)
	implicit none
		double precision, dimension(P+1,N+1,Q) :: m, dm
		m(:,1,:) = -N*dm(:,1,:)
		m(:,2:N,:) = N*(dm(:,1:N-1,:) - dm(:,2:N,:))
		m(:,N+1,:) = N*dm(:,N,:)		
	end function dxS

!! Adjoint de la dérivation en y
	function dyS(dm) result(m)
	implicit none
		double precision, dimension(P+1,N+1,Q) :: m, dm
		m(1,:,:) = -P*dm(1,:,:)
		m(2:P,:,:) = P*(dm(1:P-1,:,:) - dm(2:P,:,:))
		m(P+1,:,:) = P*dm(P,:,:)		
	end function dyS

!! Adjoint de la dérivation en t 
	function dtS(df) result(f)
	implicit none
		double precision, dimension(P+1,N+1,Q) :: f, df
		f(:,:,1) = -Q*df(:,:,1)
		f(:,:,2:Q-1) = Q*(df(:,:,1:Q-2) - df(:,:,2:Q-1))
		f(:,:,Q) = Q*df(:,:,Q-1)
	end function dtS		

!! Opérateur A
	function A(w) result(aw)
	implicit none
		double precision, dimension(P+1,N+1,Q,3) :: w
		double precision, dimension(P+1,N+1,Q+2) :: aw
		aw(:,:,1:Q) = -dxS(w(:,:,:,1)) - dyS(w(:,:,:,2)) + dt(w(:,:,:,3))
		aw(:,:,Q+1) = w(:,:,1,3)
		aw(:,:,Q+2) = w(:,:,Q,3)
	end function A

!! Adjoint de A
	function AS(aw) result(w)
	implicit none
		double precision, dimension(P+1,N+1,Q,3) :: w
		double precision, dimension(P+1,N+1,Q+2) :: aw		
		w(:,:,:,1) = -dx(aw(:,:,1:Q))
		w(:,:,:,2) = -dy(aw(:,:,1:Q))
		w(:,:,:,3) = dtS(aw(:,:,1:Q)) 
		w(:,:,1,3) = w(:,:,1,3) + aw(:,:,Q+1)
		w(:,:,Q,3) = w(:,:,Q,3) + aw(:,:,Q+2)
	end function AS
end program transport
