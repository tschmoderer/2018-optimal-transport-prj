module helpers
implicit none
	
contains 
!! Gaussienne 
	function gauss(mu,sigma,N) result(f)
		implicit none
		real :: mu, sigma
		integer :: N 
		double precision, dimension(1:N+1) :: f
		integer :: i
		do i = 1,N+1 
			f(i) = exp(-0.5*(((i-1)/(1.0*N) - mu)/sigma)**2)
		end do
	end function gauss 
	
!! Normalise 	
	function normalise(f,N) result(res)
		implicit none
		integer :: N
		double precision, dimension(1:N+1) :: f, res
		res = f/sum(f)
	end function normalise
	
!! Dérivation selon x
	function dx(m,N,Q) result(dm)
		implicit none
		integer :: N,Q
		double precision, dimension(1:Q+1,1:N+1) :: m, dm    
		dm = 0;
		dm(:,1:N) = N*(m(:,2:N+1) - m(:,1:N))
	end function dx

!! Dérivation selon t
	function dt(f,N,Q) result(df)
		implicit none
		integer :: N,Q
		double precision, dimension(1:Q+1,1:N+1)  :: f,df
		df = 0;
		df(1:Q,:) = Q*(f(2:Q+1,:) - f(1:Q,:))
	end function dt

!! Adjoint de dx
	function dxS(dm,N,Q) result(m)
		integer :: N,Q
		double precision, dimension(1:Q+1,1:N+1) :: dm, m
		m(:,1)   = -dm(:,1)
		m(:,2:N) = dm(:,1:N-1) - dm(:,2:N)
		m(:,N+1) = dm(:,N)
		m(:,1:N+1) = N*m(:,1:N+1) 
	end function dxS

!! Adjoint de dt 
	function dtS(df,N,Q) result(f)
		integer :: N,Q
		double precision, dimension(1:Q+1,1:N+1) :: df, f
		f(1,:)   = -df(1,:)
		f(2:Q,:) = df(1:Q-1,:) - df(2:Q,:)
		f(Q+1,:) = df(Q,:) 
		f(1:Q+1,:) = Q*f(1:Q+1,:)
	end function dtS

!! Opérateur A
	function A(w,N,Q) result(Aw)
		implicit none
		integer :: N,Q;
		double precision, dimension(1:Q+1,1:N+1,2) :: w 
		double precision, dimension(1:Q+3,1:N+1) :: Aw

		Aw(1:Q+1,:) = -dxS(w(:,:,1),N,Q) + dt(w(:,:,2),N,Q)
		Aw(Q+2,:) = w(Q+1,:,2)
		Aw(Q+3,:) = w(1,:,2)
	end function A

!! Opérateur adjoint de A
	function AS(Aw,N,Q) result(ASAw)
		implicit none
		integer :: N,Q;
		double precision, dimension(1:Q+3,1:N+1) :: Aw
		double precision, dimension(1:Q+1,1:N+1,2) :: ASAw 

		ASAw(:,:,1) = -dx(Aw(1:Q+1,:),N,Q)
		ASAw(:,:,2) = dtS(Aw(1:Q+1,:),N,Q)

		ASAw(1,:,2) = ASAw(1,:,2) + Aw(Q+3,:);
		ASAw(Q+1,:,2) = ASAw(Q+1,:,2) + Aw(Q+2,:);
	end function AS
	
!! Function aplatir
	function flat(x,N,Q) result(fx)
		implicit none
		integer :: N,Q;
		double precision, dimension(1:Q+3,1:N+1) :: x
		double precision, dimension(1:(Q+3)*(N+1)) :: fx
		fx = reshape(x,(/(Q+3)*(N+1)/)) ! une colonne pleins de lignes
	end function flat

!! Fonction reshape
	function resh(x,N,Q) result(rx)
		implicit none
		integer :: N,Q;
		double precision, dimension(1:(Q+3)*(N+1)) :: x
		double precision, dimension(1:Q+3,1:N+1):: rx
		rx = reshape(x,(/Q+3,N+1/));
	end function resh	

!! Fonction gradient conjugué
	function cg(b,N,Q) result(x)
		implicit none
		integer :: N,Q 
		double precision, dimension(1:(Q+3)*(N+1)) :: b, x, r, p, Ap 
		double precision :: alpha, rold, rnew
		integer :: i
		x = 0; 
		r = b - flat(A(AS(resh(x,N,Q),N,Q),N,Q),N,Q)
		p = r
		rold = sum(r*r)
		do i = 1,(Q+3)*(N+1)
			Ap = flat(A(AS(resh(p,N,Q),N,Q),N,Q),N,Q)
			alpha = rold/sum(p*Ap)
			x = x + alpha*p
			r = r - alpha*Ap
			rnew = sum(r*r)
			if (dsqrt(rnew) .LT. 1e-10) then 
				exit 
			end if
			p = r + (rnew/rold)*p
			rold = rnew 
		end do
	end function cg	
end module helpers




program test 
use helpers
    implicit none 
    integer, parameter :: N = 50, Q = 30;
	double precision, dimension(1:Q+1,1:N+1) :: rm, dxrm, rdxrm, dxSrdxrm
	double precision, dimension(1:Q+1,1:N+1) :: rf, dtrf, rdtrf, dtSrdtrf
	double precision, dimension(1:Q+1,1:N+1,2) :: rw, ASrArw
	double precision, dimension(1:Q+3,1:N+1) :: Arw, rArw
	
	! variables projection 
	double precision, dimension(1:N+1) :: f0, f1
	double precision, dimension(1:Q+3,1:N+1) :: y
	double precision, dimension(1:Q+1,1:N+1,2) :: w, pC
	double precision :: epsilon 
	
	call random_number(rm)
	call random_number(rdxrm)
	call random_number(rf)
	call random_number(rdtrf)
	
	dxrm     = dx(rm,N,Q)
	dxsrdxrm = dxS(rdxrm,N,Q)
	
	dtrf     = dt(rf,N,Q)
	dtSrdtrf = dtS(rdtrf,N,Q)

    print *, 'divergence m : ', sum(rm*dxSrdxrm - dxrm*rdxrm)
    print *, 'dérivée en t : ', sum(rf*dtSrdtrf - dtrf*rdtrf)

	call random_number(rw)
	call random_number(rArw)
	
	Arw    = A(rw,N,Q)
	ASrArw = AS(rArw,N,Q)
	
	print *, 'test A(AS) : ', sum(rw*ASrArw) - sum(rArw*Arw) 
  
	!! test projection
	epsilon = 1e-10
	f0 = normalise(epsilon + gauss(0.5,0.05,N),N)
	f1 = normalise(epsilon + gauss(0.5,0.05,N),N)
	y(1:Q+1,:) = 0
	y(Q+2,:) = f0
	y(Q+3,:) = f1
	
	call random_number(w)
	print *, 'error before projection : ', sum(dabs(A(w,N,Q) - y))/sum(y)	
	pC = w + AS(resh(cg(flat(y-A(w,N,Q),N,Q),N,Q),N,Q),N,Q) 
	print *, 'error after projection : ', sum(dabs(A(pC,N,Q) - y))/sum(y)
end program




