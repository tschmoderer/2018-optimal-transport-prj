program transport
    implicit none
    integer, parameter :: N = 31, P = 31, Q = 31, niter = 10000
    double precision, parameter :: eps = 1e-10, alpha = 1.98, g = 1./230, b = 0.5
    double precision, parameter :: pi = 4.D0*DATAN(1.D0)
    double precision, dimension(N+1,P+1) :: f0, f1
    double precision, dimension(N+1,P+1,Q+1,3) :: zV = 0, wV0 = 0, wV1 = 0
    double precision, dimension(N+2,P+2,Q+2,3) :: zU = 0, wU0 = 0, wU1 = 0, tmp = 0
	integer, dimension(N+1,P+1,Q+1) :: obstacle = 0
    double precision, dimension(niter) :: cout, minF, divV
    double precision :: t
    integer :: i, k, l 
  	character(10) :: charI;
	 
	f0 = normalise(eps + gauss(0.2d0,0.2d0,0.05d0))
    f1 = normalise(eps + gauss(0.8d0,0.8d0,0.05d0))! + gauss(0.8d0,0.5d0,0.05d0)) 
    
!    open(1,file='fourati.dat')
!    open(2,file='knippel.dat')
!    do i = 1,N+1
!		read(1,*) f0(i,:)
!		read(2,*) f1(i,:)
!    end do
!    close(1)
!    close(2)    
    
!    f0 = normalise(eps + f0)
!    f1 = normalise(eps + f1)

	do k = 1,Q+1
		open(1,file='maze.dat')
		do i = 1,N+1
!			read(1,*) obstacle(i,:,k)
		end do
		close(1)
    end do 
!    f0 = normalise(eps + gauss(10d0/128d0,26d0/128d0,0.05d0))
!    f1 = normalise(eps + gauss(119d0/128d0,103d0/128d0,0.05d0)) 

    
	do i = 1,Q+2
		t = (i-1)/(1.*(Q+1))
		wU0(1:P+1,1:N+1,i,3) = (1-t)*f1 + t*f0
    end do 
    wV0 = interp(wU0)
    zU = wU0; zV = wV0
    
    
    do i = 1,niter
		! A - DR
		wU1 = wU0 + alpha*(projC(2*zU - wU0) - zU)
		wV1 = wV0 + alpha*(proxJ(2*zV - wV0) - zV)
		zU  = projCs(wU1,wV1)
		zV  = interp(zU)
		! A - DR'
!		tmp = projCs(2*zU - wU0,2*zV - wV0)
!		wU1 = wU0 + alpha*(tmp- zU)
!		wV1 = wV0 + alpha*(interp(tmp)- zV)
!		zU  = projC(wU1)
!		zV  = proxJ(wV1)
		
		wU0 = wU1
		wV0 = wV1
		
		! record data
		cout(i) = J(zV)
		minF(i) = minval(zV(:,:,:,3))
		divV(i) = sum(div(Zu)**2)
        
        if (modulo(i,100) .EQ. 0) print *, i, cout(i), divV(i), sum(zV - interp(zU))
    end do 
  
    open(1,file='results/data.dat');
    write(1,*) "# ", "iter ", "energie ", "minF ", "divV "
    do i = 1,niter
		write(1,*) i, cout(i), minF(i), divV(i)
    end do
    close(1)  
    
    open(1,file='results/transport.dat');
    write(1,*) "# ", "X ", "Y ", "T ", "Z "
    do i = 1,N+1 ! y direction
        do k = 1,P+1 ! x direction 
					do l = 1,Q+2 ! t direction
           	write(1,*) i, k, l,  zU(i,k,l,3)
           end do
        end do
   	end do
    close(1)
    open(1,file='results/vitesse1.dat');
    write(1,*) "# ", "X ", "Y ", "T ", "Z "
    do i = 1,N+2 ! y direction
        do k = 1,P+1 ! x direction 
					do l = 1,Q+1 ! t direction
           	write(1,*) i, k, l,  zU(i,k,l,1)
           end do
        end do
   	end do
    close(1)    
    open(1,file='results/vitesse2.dat');
    write(1,*) "# ", "X ", "Y ", "T ", "Z "
    do i = 1,N+1 ! y direction
        do k = 1,P+2 ! x direction 
					do l = 1,Q+1 ! t direction
           	write(1,*) i, k, l,  zU(i,k,l,1)
           end do
        end do
   	end do
    close(1)    

	do l = 1,Q+2
		write(charI,'(I5.5)') l
		open(1,file='results/Transport/'//trim(charI)//'.dat'); 
		write(1,*) "# ", "X ", "Y ", "F "
		do i = 1,P+1 ! y direction
			do k = 1,N+1 ! x direction 
				write(1,*) (k-1)/(1.0*N), (i-1)/(1.0*P),  zU(i,k,l,3)
			end do
			!write(1,*)  zU(i,:,l,3)
		end do
		close(1)
	end do    
    
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

	open(8,file='results/plot.gnu'); 
	write(8,*) 'set dgrid3d ', P+1, ',', N+1
	write(8,*) 'set zr [', minval(minF) , ':', maxval(zU(:,:,:,3)), ']'
	close(8);
  
    contains

!! Gauss
    function gauss(muX,muY,sigma) result(f) 
		implicit none
        double precision :: muX, muY, sigma, x, y
        double precision, dimension(N+1,P+1) :: f
        integer :: i,j
        do i = 1,N+1
			do j = 1,Q+1
				x = (i-1)/(1.0*N)
				y = (j-1)/(1.0*Q)
				f(i,j) = exp(-0.5*((x-muX)**2 + (y-muY)**2)/sigma**2)
			end do
        end do 
    end function gauss
	
!! Normalise
    function normalise(f) result(nf) 
	implicit none
	    double precision, dimension(N+1,P+1) :: f, nf
        nf = f/sum(f)
    end function

!! Cout 
    function J(w) result(c) 
    implicit none
        double precision, dimension(N+1,P+1,Q+1,3) :: w
        double precision :: c
        c = 0.5*sum((w(:,:,:,1)**2 + w(:,:,:,2)**2)/max(w(:,:,:,3),eps,1e-10)**b)
    end function J 

!! Proximal de J 
    function proxJ(w) result(pw)
    implicit none
        double precision, dimension(N+1,P+1,Q+1,3) :: w, pw
        double precision, dimension(N+1,P+1,Q+1,2) :: mt
        double precision, dimension(N+1,P+1,Q+1)   :: ft, x0, x1, poly, dpoly
        double precision, dimension(N+1,P+1,Q+1)   :: d0, d1, theta, x2, x3, a
        integer :: k
        
        mt = w(:,:,:,1:2); ft = w(:,:,:,3);
        
        if (b .EQ. 0) then ! Interpolation L2
			x1 = ft
        else if (b .EQ. 1) then ! Transport
        	! résolution polynôme de degré 3
        	! http://www.aip.de/groups/soe/local/numres/bookfpdf/f5-6.pdf
			d0 = ((2*g - ft)**2 - 3.*(g**2 - 2*g*ft))/9.
			d1 = (2*(2*g - ft)**3 - 9*(2*g - ft)*(g**2 - 2*g*ft) + 27*(-ft*g**2 -0.5*g*(mt(:,:,:,1)**2 + mt(:,:,:,2)**2)))/54.
        	
        	where (d1**2 - d0**3 .LT. 0) 
				theta = dacos(d1/dsqrt(d0**3))
				x1 = -2*dsqrt(d0)*dcos(theta/3.) - (2*g - ft)/3.
				x2 = -2*dsqrt(d0)*dcos((theta+2*pi)/3.) - (2*g - ft)/3.
				x3 = -2*dsqrt(d0)*dcos((theta-2*pi)/3.) - (2*g - ft)/3.
				x1 = max(x1,x2,x3)
			else where
				a = -(d1/dabs(d1))*(dabs(d1) + dsqrt(d1**2 - d0**3))**(1/3.)
				where ( a .NE. 0)
					x1 = (a + d0/a) -(2*g - ft)/3.
				else where
					x1 = -(2*g - ft)/3.
				end where
        	end where
        	       	
       	else ! entre les deux
        		 x0 = 1; x1 = 2; k = 0;

		      do while (maxval(dabs(x0-x1)) .GT. 1e-5  .AND. k .LT. 50)
		          x0 = x1
					poly = x0**(1.0-b)*(x0-ft)*((x0**b+g)**2)-0.5*b*g*(mt(:,:,:,1)**2 + mt(:,:,:,2)**2)
					dpoly = (1.0-b)*x0**(-b)*(x0-ft)*((x0**b+g)**2) + x0**(1-b)*((x0**b+g)**2 +2*b*(x0-ft)*x0**(b-1)*(x0**b+g) )
							
						where (x0 .GT. eps) 
							x1 = x0 - poly/dpoly
						else where 
							x1 = eps		
						end where
		          k = k+1
		      end do
        	
        end if   

        where ((x1 .LT. eps) .OR. (obstacle .GT. 0)) x1 = eps
        
        pw(:,:,:,1) = (x1**b)*mt(:,:,:,1)/(x1**b+g) 
        pw(:,:,:,2) = (x1**b)*mt(:,:,:,2)/(x1**b+g) 
        pw(:,:,:,3) = x1
    end function proxJ
    
!! Interpolation 
	function interp(U) result(V)
	implicit none 
		double precision, dimension(N+2,P+2,Q+2,3) :: U
		double precision, dimension(N+1,P+1,Q+1,3) :: V
		V(:,:,:,1) = U(1:N+1,1:P+1,1:Q+1,1) + U(2:N+2,1:P+1,1:Q+1,1)
		V(:,:,:,2) = U(1:N+1,1:P+1,1:Q+1,2) + U(1:N+1,2:P+2,1:Q+1,2)
		V(:,:,:,3) = U(1:N+1,1:P+1,1:Q+1,3) + U(1:N+1,1:P+1,2:Q+2,3)
		V = 0.5*V
	end function interp

!! Interpolation adjoint 
	function interpAdj(V) result(U)
		double precision, dimension(N+2,P+2,Q+2,3) :: U
		double precision, dimension(N+1,P+1,Q+1,3) :: V
		U = 0
		
		U(1    ,1:P+1,1:Q+1,1) = V(1,:,:,1)
		U(2:N+1,1:P+1,1:Q+1,1) = V(1:N,:,:,1) + V(2:N+1,:,:,1)
		U(N+2  ,1:P+1,1:Q+1,1) = V(N+1,:,:,1)
		
		U(1:N+1,1    ,1:Q+1,2) = V(:,1,:,2)
		U(1:N+1,2:P+1,1:Q+1,2) = V(:,1:P,:,2) + V(:,2:P+1,:,2)
		U(1:N+1,P+2  ,1:Q+1,2) = V(:,P+1,:,2)
		
		U(1:N+1,1:P+1,1    ,3) = V(:,:,1,3)
		U(1:N+1,1:P+1,2:Q+1,3) = V(:,:,1:Q,3) + V(:,:,2:Q+1,3)
		U(1:N+1,1:P+1,Q+2  ,3) = V(:,:,Q+1,3)
		
		U = 0.5*U 	
	end function interpAdj

!! Projection sur Cs
	function projCs(U,V) result(pU)
	implicit none
		double precision, dimension(N+2,P+2,Q+2,3) :: U, b, r, dir, Ip, pU
		double precision, dimension(N+1,P+1,Q+1,3) :: V
		double precision :: alpha, rnew, rold
		integer :: i 
		b = U + interpAdj(V)
		pU = 0
		r = b - pU - interpAdj(interp(pU))
		dir = r
		rold = sum(r*r)
		do i = 1,3*(N+2)*(P+2)*(Q+2)
			Ip = dir + interpAdj(interp(dir))
			alpha = rold/sum(dir*Ip)
			pU = pU + alpha*dir
			r = r - alpha*Ip
			rnew = sum(r*r) 
			if (dsqrt(rnew) .LT. 1e-10) exit
			dir = r + (rnew/rold)*dir
			rold = rnew
		end do 
	end function projCs

!! Divergence 
	function div(U) result(D) 
	implicit none
		double precision, dimension(N+2,P+2,Q+2,3) :: U
		double precision, dimension(N+1,P+1,Q+1) :: D
		D =     (N+1)*(U(2:N+2,1:P+1,1:Q+1,1) - U(1:N+1,1:P+1,1:Q+1,1))
		D = D + (P+1)*(U(1:N+1,2:P+2,1:Q+1,2) - U(1:N+1,1:P+1,1:Q+1,2))
		D = D + (Q+1)*(U(1:N+1,1:P+1,2:Q+2,3) - U(1:N+1,1:P+1,1:Q+1,3))
	end function div

!! Adjoint de la divergence 
	function divAdj(D) result(U)
	implicit none
		double precision, dimension(N+2,P+2,Q+2,3):: U
		double precision, dimension(N+1,P+1,Q+1)  :: D
		U = 0
		
		U(1    ,1:P+1,1:Q+1,1) = - D(1,:,:)
		U(2:N+1,1:P+1,1:Q+1,1) =   D(1:N,:,:) - D(2:N+1,:,:)
		U(N+2  ,1:P+1,1:Q+1,1) =   D(N+1,:,:)
		U(:,:,:,1) = (N+1)*U(:,:,:,1)
		
		U(1:N+1,1    ,1:Q+1,2) = - D(:,1,:)
		U(1:N+1,2:P+1,1:Q+1,2) =   D(:,1:P,:) - D(:,2:P+1,:)
		U(1:N+1,P+2  ,1:Q+1,2) =   D(:,P+1,:)
		U(:,:,:,2) = (P+1)*U(:,:,:,2)
		
		U(1:N+1,1:P+1,1    ,3) = - D(:,:,1)
		U(1:N+1,1:P+1,2:Q+1,3) =   D(:,:,1:Q) - D(:,:,2:Q+1)
		U(1:N+1,1:P+1,Q+2  ,3) =   D(:,:,Q+1)
		U(:,:,:,3) = (Q+1)*U(:,:,:,3)
	end function divAdj

!! Projection sur C
	function projC(U) result(pU) 
	implicit none 
		double precision, dimension(N+2,P+2,Q+2,3) :: U, pU, gf
		double precision, dimension(N+1,P+1,Q+1)   :: D, f	

		U(1:N+1,1:P+1,1  ,3) = f0
		U(1:N+1,1:P+1,Q+2,3) = f1
		
		U(1  ,1:P+1,1:Q+1,1) = 0
		U(N+2,1:P+1,1:Q+1,1) = 0
		
		U(1:N+1,1  ,1:Q+1,2) = 0
		U(1:N+1,P+2,1:Q+1,2) = 0
		
		D = div(U)
		f = poisson(-D)
		
		gf = divAdj(f)
		
		pU = U
		pU(2:N+1,1:P+1,1:Q+1,1) = pU(2:N+1,1:P+1,1:Q+1,1) + gf(2:N+1,1:P+1,1:Q+1,1)
		pU(1:N+1,2:P+1,1:Q+1,2) = pU(1:N+1,2:P+1,1:Q+1,2) + gf(1:N+1,2:P+1,1:Q+1,2)
		pU(1:N+1,1:P+1,2:Q+1,3) = pU(1:N+1,1:P+1,2:Q+1,3) + gf(1:N+1,1:P+1,2:Q+1,3)
	end function projC

!! poisson 
	function poisson(f) result(poi)
	implicit none 
		double precision, dimension(N+1,P+1,Q+1) :: f, poi, denom, fhat, uhat
		double precision, dimension(N+1) :: dn, depn
		double precision, dimension(P+1) :: dp, depp
		double precision, dimension(Q+1) :: dq, depq
		integer :: i, j
		
		do i = 1,N+1
			dn(i) = (i-1)/(1.*(N+1))
		end do
		depn = (2*dcos(pi*dn) - 2)*(N+1)**2
		do i = 1,P+1
			dp(i) = (i-1)/(1.*(P+1))
		end do
		depp = (2*dcos(pi*dp) - 2)*(P+1)**2	
		do i = 1,Q+1
			dq(i) = (i-1)/(1.*(Q+1))
		end do
		depq = (2*dcos(pi*dq) - 2)*(Q+1)**2		
		
		
		do i = 1,N+1
			do j = 1,P+1
				denom(i,j,:) = depn(i) + depp(j) + depq
			end do 
		end do
	
		where (denom .EQ. 0) denom = 1.
		
		fhat = dct3(f,N+1,P+1,Q+1)
		uhat = -fhat/denom
		poi  = idct3(uhat,N+1,P+1,Q+1)	
	end function poisson
	
	function dct(f,M) result(df)
		implicit none
		integer :: M
		double precision, dimension(M) :: f, df
		double precision, dimension(M,M) :: ADCT
		double precision, dimension(M) :: a1
		integer :: u,x
			
		a1 = dsqrt(2d0/(1.*M)); a1(1) = 1./dsqrt(1d0*M)
		do u = 1,M
			do x = 1,M
				ADCT(u,x) = a1(u)*dcos(pi*(2*x-1)*(u-1)/(2.*M))
			end do 
		end do 
		df = matmul(ADCT,f)
	end function dct
	
	function idct(df,M) result(f)
		implicit none
		integer :: M
		double precision, dimension(M) :: f, df
		double precision, dimension(M,M) :: ADCT
		double precision, dimension(M) :: a1
		integer :: u,x
			
		a1 = dsqrt(2d0/(1.*M)); a1(1) = 1./dsqrt(1d0*M)
		do u = 1,M
			do x = 1,M
				ADCT(u,x) = a1(u)*dcos(pi*(2*x-1)*(u-1)/(2.*M))
			end do 
		end do 
		f = matmul(transpose(ADCT),df)	
	end function idct	
		
	
	function dct2(f,S1,S2) result(df)
	implicit none
		integer :: S1, S2
		double precision, dimension(S1,S2) :: f, df
		double precision, dimension(S1,S1) :: ADCT
		double precision, dimension(S2,S2) :: ADCT2
		double precision, dimension(S1) :: a1
		double precision, dimension(S2) :: a2
		integer :: u,x,v,y
			
		a1 = dsqrt(2d0/(1.*S1)); a1(1) = 1./dsqrt(1d0*S1)
		a2 = dsqrt(2d0/(1.*S2)); a2(1) = 1./dsqrt(1d0*S2)
		do u = 1,S1
			do x = 1,S1
				ADCT(u,x) = a1(u)*dcos(pi*(2*x-1)*(u-1)/(2.*S1))
			end do 
		end do 
		do v = 1,S2
			do y = 1,S2
				ADCT2(v,y) = a2(v)*dcos(pi*(2*y-1)*(v-1)/(2.*S2))
			end do 
		end do 
		
		df = matmul(ADCT,matmul(f,transpose(ADCT2)))
	end function dct2
	
	function idct2(df,S1,S2) result(f)
	implicit none
		integer :: S1, S2
		double precision, dimension(S1,S2) :: f, df
		double precision, dimension(S1,S1) :: ADCT
		double precision, dimension(S2,S2) :: ADCT2
		double precision, dimension(S1) :: a1
		double precision, dimension(S2) :: a2
		integer :: u,x,v,y
			
		a1 = dsqrt(2d0/(1.*S1)); a1(1) = 1./dsqrt(1d0*S1)
		a2 = dsqrt(2d0/(1.*S2)); a2(1) = 1./dsqrt(1d0*S2)
		do u = 1,S1
			do x = 1,S1
				ADCT(u,x) = a1(u)*dcos(pi*(2*x-1)*(u-1)/(2.*S1))
			end do 
		end do 
		do v = 1,S2
			do y = 1,S2
				ADCT2(v,y) = a2(v)*dcos(pi*(2*y-1)*(v-1)/(2.*S2))
			end do 
		end do 
		
		f = matmul(transpose(ADCT),matmul(df,ADCT2))
	end function idct2
	
	function dct3(f,S1,S2,S3) result(df)
	implicit none
		integer :: S1, S2, S3
		double precision, dimension(S1,S2,S3) :: f, df
		double precision, dimension(S1,S1) :: ADCT
		double precision, dimension(S2,S2) :: ADCT2
		double precision, dimension(S3,S3) :: ADCT3
		double precision, dimension(S1) :: a1
		double precision, dimension(S2) :: a2
		double precision, dimension(S3) :: a3
		double precision, dimension(S1,S2,S3) :: tmp1, tmp2
		integer :: u,x,v,y,w,z
		integer :: i1,i2,i3
			
		a1 = dsqrt(2d0/(1.*S1)); a1(1) = 1./dsqrt(1d0*S1)
		a2 = dsqrt(2d0/(1.*S2)); a2(1) = 1./dsqrt(1d0*S2)
		a3 = dsqrt(2d0/(1.*S3)); a3(1) = 1./dsqrt(1d0*S3)
		do u = 1,S1
			do x = 1,S1
				ADCT(u,x) = a1(u)*dcos(pi*(2*x-1)*(u-1)/(2.*S1))
			end do 
		end do 
		do v = 1,S2
			do y = 1,S2
				ADCT2(v,y) = a2(v)*dcos(pi*(2*y-1)*(v-1)/(2.*S2))
			end do 
		end do 
		do w = 1,S3
			do z = 1,S3
				ADCT3(w,z) = a3(w)*dcos(pi*(2*z-1)*(w-1)/(2.*S3))
			end do 
		end do 
		
		! F (x)_1 A1
		do i2 = 1,S2
			do i3 = 1,S3
				tmp1(:,i2,i3) = matmul(ADCT,f(:,i2,i3))
			end do
		end do
		! tmp1 (x)_2 A2
		do i3 = 1,S3
			do i1 = 1,S1
				tmp2(i1,:,i3) = matmul(ADCT2,tmp1(i1,:,i3))
			end do
		end do
		! tmp2 (x)_3 A3
		do i1 = 1,S1
			do i2 = 1,S2
				df(i1,i2,:) = matmul(ADCT3,tmp2(i1,i2,:))
			end do
		end do
	end function dct3

	function idct3(df,S1,S2,S3) result(f)
	implicit none
		integer :: S1, S2, S3
		double precision, dimension(S1,S2,S3) :: f, df
		double precision, dimension(S1,S1) :: ADCT
		double precision, dimension(S2,S2) :: ADCT2
		double precision, dimension(S3,S3) :: ADCT3
		double precision, dimension(S1) :: a1
		double precision, dimension(S2) :: a2
		double precision, dimension(S3) :: a3
		double precision, dimension(S1,S2,S3) :: tmp1, tmp2
		integer :: u,x,v,y,w,z
		integer :: i1,i2,i3

		a1 = dsqrt(2d0/(1.*S1)); a1(1) = 1./dsqrt(1d0*S1)
		a2 = dsqrt(2d0/(1.*S2)); a2(1) = 1./dsqrt(1d0*S2)
		a3 = dsqrt(2d0/(1.*S3)); a3(1) = 1./dsqrt(1d0*S3)
		
		do u = 1,S1
			do x = 1,S1
				ADCT(u,x) = a1(u)*dcos(pi*(2*x-1)*(u-1)/(2.*S1))
			end do 
		end do 
		do v = 1,S2
			do y = 1,S2
				ADCT2(v,y) = a2(v)*dcos(pi*(2*y-1)*(v-1)/(2.*S2))
			end do 
		end do 
		do w = 1,S3
			do z = 1,S3
				ADCT3(w,z) = a3(w)*dcos(pi*(2*z-1)*(w-1)/(2.*S3))
			end do 
		end do 
		
		! F (x)_1 A1
		do i2 = 1,S2
			do i3 = 1,S3
				tmp1(:,i2,i3) = matmul(transpose(ADCT),df(:,i2,i3))
			end do
		end do
		! tmp1 (x)_2 A2
		do i3 = 1,S3
			do i1 = 1,S1
				tmp2(i1,:,i3) = matmul(transpose(ADCT2),tmp1(i1,:,i3))
			end do
		end do
		! tmp2 (x)_3 A3
		do i1 = 1,S1
			do i2 = 1,S2
				f(i1,i2,:) = matmul(transpose(ADCT3),tmp2(i1,i2,:))
			end do
		end do
	end function idct3
end program transport
