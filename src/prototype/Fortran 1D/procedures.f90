module procedures
implicit none
	integer, parameter :: N = 2, Q = 2
	integer, parameter :: niter = 59
    real, parameter :: eps = 1e-10, alpha = 1.0, g = 1.0
    
contains 
!! Gaussienne 
	function gauss(mu,sigma) result(f)
		implicit none
		real :: mu, sigma
		real, dimension(1:N+1) :: f
		integer :: i
		do i = 1,N+1 
			f(i) = exp(-0.5*(((i-1)/(1.0*N) - mu)/sigma)**2)
		end do
	end function gauss 
	
!! Normalise 	
	function normalise(f) result(res)
		implicit none
		real, dimension(1:N+1) :: f, res
		res = f/sum(f)
	end function normalise
	
!! Dérivation selon x
	function dx(m) result(dm)
		implicit none
		real, dimension(1:Q+1,1:N+1) :: m, dm    
		dm = 0;
		dm(:,1:N) = N*(m(:,2:N+1) - m(:,1:N))
	end function dx

!! Dérivation selon t
	function dt(f) result(df)
		implicit none
		real, dimension(1:Q+1,1:N+1)  :: f,df
		df = 0;
		df(1:Q,:) = Q*(f(2:Q+1,:) - f(1:Q,:))
	end function dt

!! Adjoint de dx
	function dxS(dm) result(m)
		real, dimension(1:Q+1,1:N+1) :: dm, m
		m(:,1)   = -dm(:,1)
		m(:,2:N) = dm(:,1:N-1) - dm(:,2:N)
		m(:,N+1) = dm(:,N)
		m(:,1:N+1) = N*m(:,1:N+1) 
	end function dxS

!! Adjoint de dt 
	function dtS(df) result(f)
		real, dimension(1:Q+1,1:N+1) :: df, f
		f(1,:)   = -df(1,:)
		f(2:Q,:) = df(1:Q-1,:) - df(2:Q,:)
		f(Q+1,:) = df(Q,:) 
		f(1:Q+1,:) = Q*f(1:Q+1,:)
	end function dtS

!! Opérateur A
	function A(w) result(Aw)
		implicit none
		real, dimension(1:Q+1,1:N+1,2) :: w 
		real, dimension(1:Q+3,1:N+1) :: Aw

		Aw(1:Q+1,:) = -dxS(w(:,:,1)) + dt(w(:,:,2))
		Aw(Q+2,:) = w(Q+1,:,2)
		Aw(Q+3,:) = w(1,:,2)
	end function A

!! Opérateur adjoint de A
	function AS(Aw) result(ASAw)
		implicit none
		real, dimension(1:Q+3,1:N+1) :: Aw
		real, dimension(1:Q+1,1:N+1,2) :: ASAw 

		ASAw(:,:,1) = -dx(Aw(1:Q+1,:))
		ASAw(:,:,2) = dtS(Aw(1:Q+1,:))

		ASAw(1,:,2) = ASAw(1,:,2) + Aw(Q+3,:);
		ASAw(Q+1,:,2) = ASAw(Q+1,:,2) + Aw(Q+2,:);
	end function AS
	
!! Function aplatir
	function flat(x) result(fx)
		implicit none
		real, dimension(1:Q+3,1:N+1) :: x
		real, dimension(1:(Q+3)*(N+1)) :: fx
		fx = reshape(x,(/(Q+3)*(N+1)/)) ! une colonne pleins de lignes
	end function flat

!! Fonction reshape
	function resh(x) result(rx)
		implicit none
		real, dimension(1:(Q+3)*(N+1)) :: x
		real, dimension(1:Q+3,1:N+1):: rx
		rx = reshape(x,(/Q+3,N+1/));
	end function resh	

!! Fonction gradient conjugué
	function cg(b) result(x)
		implicit none
		real, dimension(1:(Q+3)*(N+1)) :: b, x, r, p, Ap 
		real :: alpha, rold, rnew
		integer :: i
		
		x = 0
		r = b - flat(A(AS(resh(x))))
		p = r
		rold = sum(r*r)
		do i = 1,(Q+3)*(N+1)
			Ap = flat(A(AS(resh(p))))
			alpha = rold/sum(p*Ap)
			x = x + alpha*p
			r = r - alpha*Ap
			rnew = sum(r*r)
			if (sqrt(rnew) .LT. 1e-10) then
				exit 
			end if
			p = r + (rnew/rold)*p
			rold = rnew 
		end do
	end function cg	

!! Projection sur C
	subroutine projC(pC,div,w,f0,f1) 
		real, dimension(1:Q+1,1:N+1,2), intent(in) :: w 
		real, dimension(1:Q+1,1:N+1,2), intent(out) :: pC
		real, dimension(1:N+1), intent(in) :: f0, f1
		real, intent(out) :: div
		real, dimension(1:Q+3,1:N+1) :: y

		y(1:Q+1,:) = 0; y(Q+2,:) = f0; y(Q+3,:) = f1;
		
		pC = w + AS(resh(cg(flat(y-A(w))))) 
		div = sum((A(pC)-y)**2)/sum(y**2)
	end subroutine projC	

!! le cout
	function cost(w) result(c)
		real, dimension(1:Q+1,1:N+1,2) :: w
		real:: c
		c = sum(w(:,:,1)**2/w(:,:,2))
	end function cost

!! Prox de J
	function proxJ(w) result(Pw)
		implicit none
		real, dimension(1:Q+1,1:N+1,2) :: w, Pw
		real, dimension(1:Q+1,1:N+1) :: mt, ft, x0, x1, poly, dpoly
		integer :: i, j, k = 0
		
		x0 = 1000; x1 = 2000; 
		mt = w(:,:,1); ft = w(:,:,2);

		do while ((sum(abs(x0-x1)) .GT. 1e-10) .AND. (k .LT. 1500))
			x0 = x1
			poly = (x0-ft)*(x0+g)*(x0+g)-0.5*g*mt*mt
			dpoly = 2*(x0+g)*(x0-ft) + (x0+g)*(x0+g)
			x1 = x0 - poly/dpoly
			k = k+1
		end do 
		Pw(:,:,2) = x1
		Pw(:,:,1) = x1*mt/(x1+g)
		
		do i = 1,Q+1
			do j = 1,N+1
				if (x1(i,j) .LT. 0) then
					Pw(i,j,1) = 0
					Pw(i,j,2) = 0
				end if
			end do
		end do
	end function proxJ
end module procedures
