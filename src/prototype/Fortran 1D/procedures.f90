module procedures
implicit none
	
contains 
!! Gaussienne 
	function gauss(mu,sigma,N) result(f)
		implicit none
		double precision :: mu, sigma
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
		
		x = 0
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

!! Projection sur C
	subroutine projC(pC,div,w,f0,f1,N,Q) 
		integer, intent(in) :: N,Q
		double precision, dimension(1:Q+1,1:N+1,2) :: w, pC
		double precision, dimension(1:N+1), intent(in) :: f0, f1
		double precision, intent(out) :: div
		double precision, dimension(1:Q+3,1:N+1) :: y

		y(1:Q+1,:) = 0; y(Q+2,:) = f0; y(Q+3,:) = f1;
		
		pC = w + AS(resh(cg(flat(y-A(w,N,Q),N,Q),N,Q),N,Q),N,Q) 
		div = sum(dabs(A(pC,N,Q)-y))/sum(y)
	end subroutine projC	

!! le cout
	function cost(w,N,Q) result(c)
		integer :: N,Q
		double precision, dimension(1:Q+1,1:N+1,2) :: w
		double precision:: c
		c = sum(w(:,:,1)**2/w(:,:,2))
	end function cost

!! Prox de J
	function proxJ(w,g,N,Q) result(Pw)
	implicit none
		integer :: N,Q
		double precision :: g
		double precision, dimension(1:Q+1,1:N+1,2) :: w, Pw
		double precision, dimension(1:Q+1,1:N+1) :: mt, ft
		double precision, dimension(1:Q+1,1:N+1) :: x0, x1, poly, dpoly, ddpoly

		integer :: k = 0, i, j
		x0 = 1000; x1 = 2000
		mt = w(:,:,1); ft = w(:,:,2)
		! Newton
		do while (maxval(dabs(x0-x1)) > 1e-10 .AND. k < 1500)
			x0 = x1;
			poly  = (x0 - ft)*((x0 + g)**2)-0.5*g*(mt**2)
			dpoly = 2.0*(x0 + g)*(x0 - ft) + (x0 + g)**2
			ddpoly = 2.0*(2.0*x0-ft+g) + 2.0*(x0+g)
		!	x1 = x0 - poly/dpoly
			x1 = x0 - 2*poly*dpoly/(2*dpoly**2-poly*ddpoly)
			k = k+1;
		end do 
	!	print *, ' k = ', k 


		    open(1,file='results/x1.dat');
			write(1,*) "# ", "X ", "T ", "Z "
			do i = 1,Q+1
				do j = 1,N+1
					write(1,*) (i-1)/(1.0*Q), (j-1)/(1.0*N), x1(i,j)
				end do
			end do
			close(1)

        Pw(:,:,2) = x1;
        Pw(:,:,1) = x1*mt/(x1+g)
		do i = 1,Q+1
			do j = 1,N+1 
				if (x1(i,j) .GT. 0) then
				!	Pw(i,j,2) = x1(i,j)
				!	Pw(i,j,1) = Pw(i,j,2)*mt(i,j)/(Pw(i,j,2) + g)
				else 
				!	Pw(i,j,2) = 0
				!	Pw(i,j,1) = 1e-10
				end if
                if (x1(i,j) .LE. 0) then
                    Pw(i,j,1) = 0
                    Pw(i,j,2) = 0
                end if
			end do 
		end do


		    open(1,file='results/PW.dat');
			write(1,*) "# ", "X ", "T ", "Z "
			do i = 1,Q+1
				do j = 1,N+1
					write(1,*) (i-1)/(1.0*Q), (j-1)/(1.0*N), PW(i,j,2)
				end do
			end do
			close(1)


	end function proxJ

end module procedures
