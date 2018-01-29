program test 
use procedures
    implicit none 
    integer, parameter :: N1 = 2, Q1 = 2;
    double precision, dimension(1:Q1+1,1:N1+1,2) :: w1, pJw, pCw;
    double precision, dimension(1:N1+1) :: f01, f11;
    integer :: i,j

    integer, parameter :: N = 100, Q = 30;
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
	f0 = normalise(epsilon + gauss(0.5d0,0.05d0,N),N)
	f1 = normalise(epsilon + gauss(0.5d0,0.05d0,N),N)
	y(1:Q+1,:) = 0
	y(Q+2,:) = f0
	y(Q+3,:) = f1
	
	call random_number(w)
	print *, 'error before projection : ', sum((A(w,N,Q) - y)**2)/sum(y**2)	
	pC = w + AS(resh(cg(flat(y-A(w,N,Q),N,Q),N,Q),N,Q),N,Q) 
	print *, 'error after projection : ', sum((A(pC,N,Q) - y)**2)/sum(y**2)


    print *,
    print *, '******************************'
    print *, 'ANOTHER PART OF THE TEST'
    print *, '******************************'
    print *,

    w1(1,1,1) = 1;
    w1(1,2,1) = 2;
    w1(1,3,1) = 3;
    w1(2,1,1) = 4;
    w1(2,2,1) = 5;
    w1(2,3,1) = 6;    
    w1(3,1,1) = 7;
    w1(3,2,1) = 8;
    w1(3,3,1) = 9;

    w1(1,1,2) = 10;
    w1(1,2,2) = 11;
    w1(1,3,2) = 12;
    w1(2,1,2) = 13;
    w1(2,2,2) = 14;
    w1(2,3,2) = 15;    
    w1(3,1,2) = 16;
    w1(3,2,2) = 17;
    w1(3,3,2) = 18;
    
    pJw = proxJ(w1,1.d0,N1,Q1)  

    do j = 1,2
        do i = 1,Q+1
!            print *, pJw(i,:,j), 'ENDL'
        end do
!        print *, 'END BLOCK'
    end do 

    f01 = 1; f11 = 1;
    ! pCw = projC(w1,f01,f11,N1,Q1)
!    print *, 'Projection sur C : '
    do j = 1,2
        do i = 1,Q+1
!            print *, pCw(i,:,j), 'ENDL'
        end do
!        print *, 'END BLOCK'
    end do 
end program




