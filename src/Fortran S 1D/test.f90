program test
	implicit none 
	    integer, parameter :: N = 9, Q = 12, niter = 200
		double precision, parameter :: eps = 1e-10, alpha = 1.0, g = 1.0, b = 1
		double precision, dimension(N+1) :: f0, f1
		integer, dimension(Q+1,N+1) :: obstacle = 0


		double precision, dimension(Q+2,N+2,2) :: U, IsV, dSu
		double precision, dimension(Q+1,N+1,2) :: V, Iu
		double precision, dimension(Q+1,N+1) :: dU, d
		
		call random_number(U)
		call random_number(V)
		
		Iu = interp(U) 
		IsV = interpAdj(V)
	
		print *, ' Interpolation :', sum(U*IsV) - sum(V*Iu)
		
		call random_number(d)
		
		du = div(U)
		dSu = divAdj(d)
		print *, ' Divergence :', sum(U*dSu) - sum(d*dU)
    contains
        
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
		U(1:Q+1,1,1)     = V(:,1,1) 
		U(1:Q+1,2:N+1,1) = V(:,2:N+1,1) + V(:,1:N,1)
		U(1:Q+1,N+2,1)   = V(:,N+1,1)
		
		U(1,1:N+1,2)     = V(1,:,2)
		U(2:Q+1,1:N+1,2) = V(2:Q+1,:,2) + V(1:Q,:,2)
		U(Q+2,1:N+1,2)   = V(Q+1,:,2)
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
		U(:,:,1)         = N*U(:,:,1)
		
		U(1,1:N+1,2)     = -D(1,:) 
		U(2:Q+1,1:N+1,2) = D(1:Q,:) - D(2:Q+1,:)
		U(Q+2,1:N+1,2)   = D(Q+1,:) 
		U(:,:,2)         = Q*U(:,:,2)
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
		Au(Q+2,1:N+1)   = U(1,1:N+1,2)   ! frontieres de fbar
		Au(Q+3,1:N+1)   = U(Q+2,1:N+1,2)
	end function A 

!! Adjoint de A 
	function AS(R) result(U) 
	implicit none 
		double precision, dimension(Q+2,N+2,2) :: U
		double precision, dimension(Q+3,N+3) :: R
		U = 0
		U(1:Q+2,1:N+2,:) = divAdj(R(1:Q+1,1:N+1)) 
		U(1:Q+1,1,1) = U(1:Q+1,1,1) + R(1:Q+1,N+2)
		U(1:Q+1,N+2,1) = U(1:Q+1,N+2,1) + R(1:Q+1,N+3)
		U(1,1:N+1,2) = U(1,1:N+1,2) + R(Q+2,1:N+1)
		U(Q+2,1:N+1,2) = U(Q+2,1:N+1,2) + R(Q+3,1:N+1)
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
end program test
