program transport
    implicit none
    integer, parameter :: N = 15, P = 15, Q = 15, niter = 200
    double precision, parameter :: eps = 1e-10, alpha = 1.0, g = 1.0, b = 1.0
    double precision, dimension(P+1,N+1) :: f0, f1
    double precision, dimension(P+1,N+1,Q+1,3) :: zV = 0, wV0 = 0, wV1 = 0
    double precision, dimension(P+2,N+2,Q+2,3) :: zU = 0, wU0 = 0, wU1 = 0
	integer, dimension(P+1,N+1,Q+1) :: obstacle = 0
    double precision, dimension(niter) :: cout, minF
    integer :: i, k, l 
  	character(10) :: charI;
  	
  	!! var pour les tests adj
!  	double precision, dimension(P+2,N+2,Q+2,3) :: U, pU, ISV, dSD, ASrA
!  	double precision, dimension(P+1,N+1,Q,3) :: V, IU
!  	double precision, dimension(P+1,N+1,Q) :: D, dU
!   double precision, dimension(P+3,N+3,Q+2) :: Au, rA
    
    f0 = normalise(eps + gauss(0.2d0,0.2d0,0.05d0))
    f1 = normalise(eps + gauss(0.8d0,0.8d0,0.05d0))

    !! test adjoint 
!    call random_number(U)
!    IU  = interp(U)
!    call random_number(V)
!    ISV = interpAdj(V)
    
!    print *, 'Interp adjoint     ', sum(U*ISV), sum(V*IU) 
    
!    call random_number(U) 
!    dU = div(U)
!    call random_number(D)
!    dSD = divAdj(D)
    
!    print *, 'Divergence adjoint ', sum(U*dSD), sum(D*dU) 

!	call random_number(U) 
!	Au = A(U)
!	call random_number(rA)
!	ASrA = AS(rA)
	
!	print *, 'Opérateur A adjoint', sum(U*ASrA), sum(Au*rA) 

!	call random_number(U)
!	pU = projC(U)
!	print *, 'avant projection : ', sum(div(U))
!	print *, 'après projection : ', sum(div(pU))

!	call exit()

do i = 1,niter
		wU1 = wU0 + alpha*(projC(2*zU - wU0) - zU)
		wV1 = wV0 + alpha*(proxJ(2*zV - wV0) - zV)
		zU  = projCs(wU1,wV1)
		zV  = interp(zU)

		wU0 = wU1
		wV0 = wV1
		
        cout(i) = J(zV)
        minF(i) = minval(zV(:,:,:,3))

        if (modulo(i,1) .EQ. 0) print *, i, cout(i)
    end do 
    
	open(1,file='results/transport.dat');
    write(1,*) "# ", "X ", "Y ", "T ", "Z "
    do i = 1,P+1 ! y direction
        do k = 1,N+1 ! x direction 
			do l = 1,Q+1 ! t direction
           write(1,*) (k-1)/(1.0*N), (i-1)/(1.0*p), (l-1)/(1.0*Q),  zV(i,k,l,3)
           end do
        end do
    end do
    close(1)

	do l = 1,Q+1 
		write(charI,'(I5.5)') l
		open(1,file='results/Transport/'//trim(charI)//'.dat'); 
		write(1,*) "# ", "X ", "Y ", "T ", "Z "
		do i = 1,P+1 ! y direction
			do k = 1,N+1 ! x direction 
				write(1,*) (k-1)/(1.0*N), (i-1)/(1.0*P),  zV(i,k,l,3)
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
    write(1,*) "# ", "X ", "Y "
    do i = 1,P+1
		do k = 1,N+1
			write(1,*) (k-1)/(1.0*N), (i-1)/(1.0*N), f0(i,k)
        end do
    end do
    close(1)

    open(1,file='results/f1.dat');
    write(1,*) "# ", "X ", "Y "
    do i = 1,P+1
		do k = 1,N+1
			write(1,*) (k-1)/(1.0*N), (i-1)/(1.0*N), f1(i,k)
        end do
    end do
    close(1)

	open(8,file="results/plot.gnu"); 
	write(8,*) 'set dgrid3d ', P+1, ',', N+1
	write(8,*) 'set zr [', minval(minF) , ':', maxval(zV(:,:,:,3)), ']'
	close(8);

    contains

!! Gauss
    function gauss(muX,muY,sigma) result(f) 
		implicit none
        double precision :: muX, muY, sigma, x, y
        double precision, dimension(P+1,N+1) :: f
        integer :: i,j
        do i = 1,P+1
			do j = 1,N+1
				x = (j-1)/(1.0*N)
				y = (i-1)/(1.0*P)
				f(i,j) = exp(-0.5*((x-muX)**2 + (y-muY)**2)/sigma**2)
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
        double precision, dimension(P+1,N+1,Q+1,3) :: w
        double precision :: c
        c = 0.5*sum(sum(w(:,:,:,1:2)**2)/max(w(:,:,:,3),eps,1e-10)**b)
    end function J 

!! Proximal de J 
    function proxJ(w) result(pw)
    implicit none
        double precision, dimension(P+1,N+1,Q+1,3) :: w, pw
        double precision, dimension(P+1,N+1,Q+1,2) :: mt
        double precision, dimension(P+1,N+1,Q+1) :: ft, x0, x1, poly, dpoly
        integer :: k
        
        mt = w(:,:,:,1:2); ft = w(:,:,:,3);
        x0 = 1; x1 = 2; k = 0;

        do while (maxval(dabs(x0-x1)) .GT. 1e-10  .AND. k .LT. 1500)
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
        where (obstacle .GT. 0) x1 = eps
        pw(:,:,:,3) = x1
        pw(:,:,:,1) = (x1**b)*mt(:,:,:,1)/(x1**b+g) 
        pw(:,:,:,2) = (x1**b)*mt(:,:,:,2)/(x1**b+g) 
    end function proxJ
    
!! Interpolation 
	function interp(U) result(V)
	implicit none 
		double precision, dimension(P+2,N+2,Q+2,3) :: U
		double precision, dimension(P+1,N+1,Q+1,3) :: V
		V(:,:,:,1) = U(1:P+1,1:N+1,1:Q+1,1) + U(1:P+1,2:N+2,1:Q+1,1)
		V(:,:,:,2) = U(1:P+1,1:N+1,1:Q+1,2) + U(2:P+2,1:N+1,1:Q+1,2)
		V(:,:,:,3) = U(1:P+1,1:N+1,1:Q+1,3) + U(1:P+1,1:N+1,2:Q+2,3)
		V = 0.5*V
	end function interp

!! Interpolation adjoint 
	function interpAdj(V) result(U)
		double precision, dimension(P+2,N+2,Q+2,3) :: U
		double precision, dimension(P+1,N+1,Q+1,3) :: V
		
		U = 0
		U(1:P+1,1,1:Q+1,1)     = V(:,1,:,1)
		U(1:P+1,2:N+1,1:Q+1,1) = V(:,2:N+1,:,1) + V(:,1:N,:,1)
		U(1:P+1,N+2,1:Q+1,1)   = V(:,N+1,:,1)
		
		U(1,1:N+1,1:Q+1,2)     = V(1,:,:,2)
		U(2:P+1,1:N+1,1:Q+1,2) = V(2:P+1,:,:,2) + V(1:P,:,:,2)
		U(P+2,1:N+1,1:Q+1,2)   = V(P+1,:,:,2)
		
		U(1:P+1,1:N+1,1,3)     = V(:,:,1,3)
		U(1:P+1,1:N+1,2:Q+1,3) = V(:,:,2:Q+1,3) + V(:,:,1:Q,3)
		U(1:P+1,1:N+1,Q+2,3)   = V(:,:,Q+1,3)
				
		U = 0.5*U 	
	end function interpAdj
	
!! Projection sur Cs
	function projCs(U,V) result(pU)
	implicit none
		double precision, dimension(P+2,N+2,Q+2,3) :: U, b, r, dir, Idir, pU
		double precision, dimension(P+1,N+1,Q+1,3) :: V
		double precision :: alpha, rnew, rold
		integer :: i 
		b = U + interpAdj(V)
		pU = 0
		r = b - pU - interpAdj(interp(pU))
		dir = r
		rold = sum(r*r)
		do i = 1,3*(P+2)*(N+2)*(Q+2)
			Idir = dir + interpAdj(interp(dir))
			alpha = rold/sum(dir*Idir)
			pU = pU + alpha*dir
			r = r - alpha*Idir
			rnew = sum(r*r) 
			if (dsqrt(rnew) .LT. 1e-10) exit
			dir = r + (rnew/rold)*dir
			rold = rnew
		end do 
	end function projCs	
	
!! Divergence 
	function div(U) result(D) 
	implicit none
		double precision, dimension(P+2,N+2,Q+2,3) :: U
		double precision, dimension(P+1,N+1,Q+1) :: D
		D = (N)*(U(1:P+1,2:N+2,1:Q+1,1) - U(1:P+1,1:N+1,1:Q+1,1))
		D = D + (P)*(U(2:P+2,1:N+1,1:Q+1,2) - U(1:P+1,1:N+1,1:Q+1,2))
		D = D + (Q)*(U(1:P+1,1:N+1,2:Q+2,3) - U(1:P+1,1:N+1,1:Q+1,3))
	end function div	

!! Adjoint de la divergence 
	function divAdj(D) result(U)
	implicit none
		double precision, dimension(P+2,N+2,Q+2,3) :: U
		double precision, dimension(P+1,N+1,Q+1) :: D
		U = 0
		
		U(1:P+1,1,1:Q+1,1)     = -D(:,1,:)
		U(1:P+1,2:N+1,1:Q+1,1) = D(:,1:N,:) - D(:,2:N+1,:)
		U(1:P+1,N+2,1:Q+1,1)   = D(:,N+1,:)
		U(:,:,:,1) = (N)*U(:,:,:,1)
		
		U(1,1:N+1,1:Q+1,2)     = -D(1,:,:)
		U(2:P+1,1:N+1,1:Q+1,2) = D(1:P,:,:) - D(2:P+1,:,:)
		U(P+2,1:N+1,1:Q+1,2)   = D(P+1,:,:)
		U(:,:,:,2) = (P)*U(:,:,:,2)
		
		U(1:P+1,1:N+1,1,3)   = -D(:,:,1)
		U(1:P+1,1:N+1,2:Q+1,3) = D(:,:,1:Q) - D(:,:,2:Q+1)
		U(1:P+1,1:N+1,Q+2,3) = D(:,:,Q+1)
		U(:,:,:,3) = (Q)*U(:,:,:,3)
	end function divAdj

!! Opérateur A 
	function A(U) result(Au)
	implicit none
		double precision, dimension(P+2,N+2,Q+2,3) :: U
		double precision, dimension(P+3,N+3,Q+3) :: Au
		Au = 0
		Au(1:P+1,1:N+1,1:Q+1) = div(U)

		Au(1:P+1,N+2,1:Q+1)   = U(1:P+1,1,1:Q+1,1) ! frontiere de mbar1
		Au(1:P+1,N+3,1:Q+1)   = U(1:P+1,N+2,1:Q+1,1)
		
		Au(P+2,1:N+1,1:Q+1)   = U(1,1:N+1,1:Q+1,2) ! frontiere de mbar2
		Au(P+3,1:N+1,1:Q+1)   = U(P+2,1:N+1,1:Q+1,2)
		
		Au(1:P+1,1:N+1,Q+2)   = U(1:P+1,1:N+1,1,3) ! frontieres de fbar
		Au(1:P+1,1:N+1,Q+3)   = U(1:P+1,1:N+1,Q+2,3)
	end function A 		
	
!! Adjoint de A 
	function AS(R) result(U) 
	implicit none 
		double precision, dimension(P+2,N+2,Q+2,3) :: U
		double precision, dimension(P+3,N+3,Q+3) :: R
		U = 0
		
		U(1:P+2,1:N+2,1:Q+2,:) = divAdj(R(1:P+1,1:N+1,1:Q+1)) 
		
		U(1:P+1,1,1:Q+1,1)     = U(1:P+1,1,1:Q+1,1) + R(1:P+1,N+2,1:Q+1)
		U(1:P+1,N+2,1:Q+1,1)   = U(1:P+1,N+2,1:Q+1,1) + R(1:P+1,N+3,1:Q+1)
		
		U(1,1:N+1,1:Q+1,2)     = U(1,1:N+1,1:Q+1,2) + R(P+2,1:N+1,1:Q+1) 
		U(P+2,1:N+1,1:Q+1,2)   = U(P+2,1:N+1,1:Q+1,2) + R(P+3,1:N+1,1:Q+1) 
		
		U(1:P+1,1:N+1,1,3)     = U(1:P+1,1:N+1,1,3) + R(1:P+1,1:N+1,Q+2)
		U(1:P+1,1:N+1,Q+2,3)   = U(1:P+1,1:N+1,Q+2,3) + R(1:P+1,1:N+1,Q+3)
	end function AS

!! Projection sur C
	function projC(U) result(pU) 
	implicit none 
		double precision, dimension(P+2,N+2,Q+2,3) :: U, pU
		double precision, dimension(P+3,N+3,Q+3) :: y, x, b, r, dir, Adir
		double precision :: alpha, rnew, rold
		integer :: i 
		
		y = 0
		y(1:P+1,1:N+1,Q+2) = f0
		y(1:P+1,1:N+1,Q+3) = f1
		
		x = 0
		b = y - A(U)
		r = b - A(AS(x))
		dir = r
		rold = sum(r*r)
		do i = 1,(P+3)*(N+3)*(Q+3) ! attention ici 
			Adir = A(AS(dir))
			alpha = rold/sum(dir*Adir)
			x = x + alpha*dir
			r = r - alpha*Adir
			rnew = sum(r*r)
			if (dsqrt(rnew) .LT. 1e-10) exit
			dir = r + (rnew/rold) * dir
			rold = rnew
		end do
		
		pU = U + AS(x)
	end function projC
end program transport
