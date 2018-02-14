program transport
    implicit none
    integer, parameter :: N = 100, Q = 100, niter = 5000
    double precision, parameter :: eps = 1e-10, alpha = 1.0, g = 1.0, b = 1.0
    double precision, dimension(N+1) :: f0, f1
    double precision, dimension(Q+1,N+1,2) :: zV = 0, wV0 = 0, wV1 = 0
    double precision, dimension(Q+2,N+2,2) :: zU = 0, wU0 = 0, wU1 = 0
		integer, dimension(Q+1,N+1) :: obstacle = 0
    double precision, dimension(niter) :: cout, minF
    integer :: i, k, l 
  	character(10) :: charI;
	 
    f0 = normalise(eps + gauss(0.2d0,0.05d0))
    f1 = normalise(eps + gauss(0.8d0,0.05d0))
    
    obstacle(5,:) = 1; obstacle(5,90:95) = 0;
    obstacle(10,:) = 1; obstacle(10,5:10) = 0;
    obstacle(15,:) = 1; obstacle(15,45:50) = 0;
    obstacle(25,:) = 1; obstacle(25,75:80) = 0;
    obstacle(30,:) = 1; obstacle(30,70:75) = 0;
    obstacle(50,:) = 1; obstacle(50,5:10) = 0;
    obstacle(90,:) = 1; obstacle(50,5:10) = 0;
    
!    f0 = normalise(eps + indicatrix(0.2d0,0.3d0))
!    f1 = normalise(eps + indicatrix(0.8d0,0.9d0))
    
!    f0 = normalise(eps + gauss(0.2d0,0.05d0))
!    f1 = normalise(eps + gauss(0.8d0,0.05d0) + 0.6*gauss(0.4d0,0.05d0))
    
    do i = 1,niter
				wU1 = wU0 + alpha*(projC(2*zU - wU0) - zU)
				wV1 = wV0 + alpha*(proxJ(2*zV - wV0) - zV)
				zU  = projCs(wU1,wV1)
				zV  = interp(zU)

				wU0 = wU1
				wV0 = wV1
		
        cout(i) = J(zV)
       
        if (modulo(i,50) .EQ. 0) print *, i, cout(i)
        minF(i) = minval(zV(:,:,2))
    end do 
    
	  open(1,file='results/transport.dat');
    write(1,*) "# ", "X ", "T ", "Z "
    do i = 1,Q+1
        do k = 1,N+1
           write(1,*) (k-1)/(1.0*N), (Q - i + 1)/(1.0*Q),  zV(i,k,2)
        end do
    end do
    close(1)
    
    do l = 1,Q+1
			write(charI,'(I5.5)') Q+2 - l
			open(1,file='results/Transport/'//trim(charI)//'.dat'); 
			write(1,*) "# ", "X ", "T "
				do k = 1,N+1 ! x direction 
					write(1,*) (k-1)/(1.0*N),  zV(l,k,2)
			end do
			close(1)
	  end do 
    
    
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
	write(8,*) 'set palette gray'
	write(8,*) 'set view 0,0'
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

!! Indicatrice 
	function indicatrix(a,b) result(f)
	implicit none
        double precision :: a, b
        double precision, dimension(N+1) :: f, x
        integer :: i
        do i = 1,N+1
            x(i) = (i-1)/(1.*N)
        end do 
        f = 0
        where ((x .GT. a) .AND. (x .LT. b)) f = 1
    end function indicatrix	
	
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
    
!! Interpolation 
	function interp(U) result(V)
	implicit none 
		double precision, dimension(Q+2,N+2,2) :: U
		double precision, dimension(Q+1,N+1,2) :: V
		V(:,:,1) = U(1:Q+1,1:N+1,1) + U(1:Q+1,2:N+2,1)
		V(:,:,2) = U(1:Q+1,1:N+1,2) + U(2:Q+2,1:N+1,2)
		V = 0.5*V
	end function interp

!! Interpolation adjoint 
	function interpAdj(V) result(U)
		double precision, dimension(Q+2,N+2,2) :: U
		double precision, dimension(Q+1,N+1,2) :: V
		U = 0
		U(1:Q+1,1,1) = V(:,1,1) 
		U(1:Q+1,2:N+1,1) = V(:,2:N+1,1) + V(:,1:N,1)
		U(1:Q+1,N+2,1) = V(:,N+1,1)
		
		U(1,1:N+1,2) = V(1,:,2)
		U(2:Q+1,1:N+1,2) = V(2:Q+1,:,2) + V(1:Q,:,2)
		U(Q+2,1:N+1,2) = V(Q+1,:,2)
		U = 0.5*U 	
	end function interpAdj

!! Projection sur Cs
	function projCs(U,V) result(pU)
	implicit none
		double precision, dimension(Q+2,N+2,2) :: U, b, r, p, Ip, pU
		double precision, dimension(Q+1,N+1,2) :: V
		double precision :: alpha, rnew, rold
		integer :: i 
		b = U + interpAdj(V)
		pU = 0
		r = b - pU - interpAdj(interp(pU))
		p = r
		rold = sum(r*r)
		do i = 1,2*(Q+2)*(N+2)
			Ip = p + interpAdj(interp(p))
			alpha = rold/sum(p*Ip)
			pU = pU + alpha*p
			r = r - alpha*Ip
			rnew = sum(r*r) 
			if (dsqrt(rnew) .LT. 1e-10) exit
			p = r + (rnew/rold)*p
			rold = rnew
		end do 
	end function projCs

!! Divergence 
	function div(U) result(D) 
	implicit none
		double precision, dimension(Q+2,N+2,2) :: U
		double precision, dimension(Q+1,N+1) :: D
		D = N*(U(1:Q+1,2:N+2,1) - U(1:Q+1,1:N+1,1)) +Q*(U(2:Q+2,1:N+1,2) - U(1:Q+1,1:N+1,2))
	end function div

!! Adjoint de la divergence 
	function divAdj(D) result(U)
	implicit none
		double precision, dimension(Q+2,N+2,2) :: U
		double precision, dimension(Q+1,N+1) :: D
		U = 0
		
		U(1:Q+1,1,1)     = -D(:,1)
		U(1:Q+1,2:N+1,1) = D(:,1:N) - D(:,2:N+1)
		U(1:Q+1,N+2,1)   = D(:,N+1)
		U(:,:,1) = N*U(:,:,1)
		
		U(1,1:N+1,2)     = -D(1,:) 
		U(2:Q+1,1:N+1,2) = D(1:Q,:) - D(2:Q+1,:)
		U(Q+2,1:N+1,2)   = D(Q+1,:) 
		U(:,:,2) = Q*U(:,:,2)
	end function divAdj

!! Opérateur A 
	function A(U) result(Au)
	implicit none
		double precision, dimension(Q+2,N+2,2) :: U
		double precision, dimension(Q+3,N+3) :: Au
		Au = 0
		Au(1:Q+1,1:N+1) = div(U)
		
		Au(1:Q+1,N+2)   = U(1:Q+1,1,1)   ! frontieres de mbar
		Au(1:Q+1,N+3)   = U(1:Q+1,N+2,1) 
		
		Au(Q+2,1:N+1)   = U(Q+2,1:N+1,2)   ! frontieres de fbar
		Au(Q+3,1:N+1)   = U(1,1:N+1,2)
	end function A 

!! Adjoint de A 
	function AS(R) result(U) 
	implicit none 
		double precision, dimension(Q+2,N+2,2) :: U
		double precision, dimension(Q+3,N+3) :: R
		U = 0
		U(1:Q+2,1:N+2,:) = divAdj(R(1:Q+1,1:N+1)) 
		
		U(1:Q+1,1,1)     = U(1:Q+1,1,1) + R(1:Q+1,N+2)
		U(1:Q+1,N+2,1)   = U(1:Q+1,N+2,1) + R(1:Q+1,N+3)
		
		U(1,1:N+1,2)     = U(1,1:N+1,2) + R(Q+3,1:N+1)
		U(Q+2,1:N+1,2)   = U(Q+2,1:N+1,2) + R(Q+2,1:N+1)
	end function AS

!! Projection sur C
	function projC(U) result(pU) 
	implicit none 
		double precision, dimension(Q+2,N+2,2) :: U, pU
		double precision, dimension(Q+3,N+3) :: y, x, b, r, p, Ap
		double precision :: alpha, rnew, rold
		integer :: i 
		y = 0
		y(Q+2,1:N+1) = f0
		y(Q+3,1:N+1) = f1
		
		x = 0
		b = y - A(U)
		r = b - A(AS(x))
		p = r
		rold = sum(r*r)
		do i = 1,(Q+3)*(N+3)
			Ap = A(AS(p))
			alpha = rold/sum(p*Ap)
			x = x + alpha*p
			r = r - alpha*Ap
			rnew = sum(r*r)
			if (dsqrt(rnew) .LT. 1e-10) exit
			p = r + (rnew/rold) * p
			rold = rnew
		end do
		
		pU = U + AS(x)
	end function projC

end program transport
