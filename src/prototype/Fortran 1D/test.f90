program test 
use procedures
    implicit none 

	real, dimension(1:Q+1,1:N+1) :: rm, dxrm, rdxrm, dxSrdxrm
	real, dimension(1:Q+1,1:N+1) :: rf, dtrf, rdtrf, dtSrdtrf
	real, dimension(1:Q+1,1:N+1,2) :: rw, ASrArw
	real, dimension(1:Q+3,1:N+1) :: Arw, rArw
	
	! variables projection 
	real , dimension(1:N+1) :: f0, f1
	real , dimension(1:Q+3,1:N+1) :: y
	real , dimension(1:Q+1,1:N+1,2) :: w, pC, pW, pW2
	real, dimension(1:(Q+3)*(N+1)) :: positif
	real :: div, J1, J2, J3
	integer :: i,j
	
	call random_number(rm)
	call random_number(rdxrm)
	call random_number(rf)
	call random_number(rdtrf)
	
	dxrm     = dx(rm)
	dxsrdxrm = dxS(rdxrm)
	
	dtrf     = dt(rf)
	dtSrdtrf = dtS(rdtrf)

    print *, 'divergence m : ', sum(rm*dxSrdxrm - dxrm*rdxrm)
    print *, 'dérivée en t : ', sum(rf*dtSrdtrf - dtrf*rdtrf)

	call random_number(rw)
	call random_number(rArw)
	
	Arw    = A(rw)
	ASrArw = AS(rArw)
	
	print *, 'test A(AS) : ', sum(rw*ASrArw) - sum(rArw*Arw) 
  
	!! test projection
	f0 = normalise(eps + gauss(0.5,0.05))
	f1 = normalise(eps + gauss(0.5,0.05))
	y(1:Q+1,:) = 0
	y(Q+2,:) = f0
	y(Q+3,:) = f1
	
	call random_number(w)
	print *, 'error before projection : ', sum((A(w) - y)**2)/sum(y**2)	
	pC = w + AS(resh(cg(flat(y-A(w))))) 
	print *, 'error after projection : ', sum((A(pC) - y)**2)/sum(y**2)

	w(1,1,1) = 1;
	w(1,2,1) = 2;
	w(1,3,1) = 3;
	w(2,1,1) = 4;
	w(2,2,1) = 5;
	w(2,3,1) = 6;
	w(3,1,1) = 7;
	w(3,2,1) = 8;
	w(3,3,1) = 9;
	
	w(1,1,2) = 10;
	w(1,2,2) = 11;
	w(1,3,2) = 12;
	w(2,1,2) = 13;
	w(2,2,2) = 14;
	w(2,3,2) = 15;
	w(3,1,2) = 16;
	w(3,2,2) = 17;
	w(3,3,2) = 18;	
	
   pw = proxJ(w)
   
	do i = 1,2
		do j = 1,N+1
	!		print *, Pw(j,:,i), 'ENDL'
		end do 
	!	print *, 'END BLOCK'
	end do
	
	f0 = 1;
	f1 = 1;
	
	call projC(pC,div,w,f0,f1)
	
!	print *, ''
!	print *, '*********************************'
!	print *, ''
	do i = 1,2
		do j = 1,N+1
	!		print *, pC(j,:,i), 'ENDL'
		end do 
	!	print *, 'END BLOCK'
	end do
	
!	print *, ''
!	print *, '*********************************'
!	print *, ''
!	print *, 'div', div
	
	call random_number(positif) 
!	print *, ''
!	print *, '*********************************'
!	print *, ''
!	print *, 'test AAS > 0 :', sum(positif*flat(A(AS(resh(positif)))))	
	
	call random_number(pW2)
!	J1 = cost(pW) +  0.5*sum((w-pW)**2)
!	J2 = cost(pW2) + 0.5*sum((w-pW2)**2)
!	J3 = cost(pW+0.0001) + 0.5*sum((w-pW-0.0001)**2)
	
!	print *, ''
!	print *, '*********************************'
!	print *, ''	
!	print *, 'test prox J :'
!	print *, 'pW : ', J1, 'pW2 :', J2, 'pW3 :', J3
	
!	print *, ''
!	print *, '*********************************'
!	print *, ''		
end program




