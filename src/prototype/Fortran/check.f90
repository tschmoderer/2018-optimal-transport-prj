!! Test des subroutines de la méthode !!

subroutine check(alpha,beta,gamma,N,Q)
    implicit none 
    double precision, intent(in) :: alpha,beta,gamma;
    integer, intent(in) :: N,Q;
    integer :: i,j;

    double precision, dimension(Q+1,N+1) :: m,f;
    double precision, dimension(Q+1,N+2) :: mbar;
    double precision, dimension(Q+2,N+1) :: fbar;

    !! Le cout de la solution !!
    double precision R;
    !! Divergence !! 
    double precision, dimension(Q+1,N+1) :: d;
    !! Boundary !!
    double precision, dimension(Q+1) :: mleft,mright;
    double precision, dimension(N+1) :: fup,fdown;
    !! Interpolsation !!
    double precision, dimension(Q+1,N+1) :: Interpm,Interpf;
    !! Prox J !!
    double precision, dimension(Q+1,N+1) :: Proxm,Proxf; 
    !! Proj C !!
    double precision, dimension(Q+1,N+2) :: projCmbar;
    double precision, dimension(Q+2,N+1) :: projCfbar;  
    !! Prox G1 !!
    double precision, dimension(Q+1,N+1) :: Vm,Vf;
    double precision, dimension(Q+1,N+2) :: Umbar;
    double precision, dimension(Q+2,N+1) :: Ufbar;
    !! Prox G2 !!
    double precision, dimension(Q+1,N+1) :: mt,ft;
    double precision, dimension(Q+1,N+2) :: mbart;
    double precision, dimension(Q+2,N+1) :: fbart;

    !!! Zone de tests !!!
    call RANDOM_NUMBER(m);
    call RANDOM_NUMBER(f);
    call RANDOM_NUMBER(mbar);
    call RANDOM_NUMBER(fbar);

!    open(0,file="../m"); write(0,*) m; close(0);
!    open(0,file="../f"); write(0,*) f; close(0);
!    open(0,file="../mbar"); write(0,*) mbar; close(0);
!    open(0,file="../fbar"); write(0,*) fbar; close(0);

!    print *, "mbar : "
!    do i = 1,Q+1
!        print *, mbar(:,i), "ENDL";
!    end do
!    print *, "fbar : "
!    do i = 1,Q+2
!        print *, fbar(:,i), "ENDL";
!    end do

    !! Check cost !!
    call cost(R,m,f,N,Q);
!    print *, "Cout de la solution : ", R;
 
    !! Check divergence !!
    call divergence(mbar,fbar,d,N,Q);
!    print *, "Divergence : ";
!    do i = 1,Q+1
!      print *, d(:,i), "ENDL";
!    end do

    !! Check Boundary !!
    call boundary(mbar,fbar,mleft,mright,fup,fdown,N,Q)
!    print *, "mleft : ", mleft; 
!    print *, "mright : ", mright;
!    print *, "fup : ", fup;
!    print *, "fdown : ", fdown;

    !! Check Interpolation 
    call interpolation(mbar,fbar,Interpm,Interpf,N,Q)
!    print *, "Interpolation m : ";
!    do i = 1,Q+1
!        print *, Interpm(:,i), "ENDL";
!    end do
!    print *, "Interpolation f : ";
!    do i = 1,Q+1
!        print *, Interpf(:,i), "ENDL";
!    end do

    !! check adjoint
    call check_adjoint(N,Q)

    !! Check Proximity operator of J
!    open(0,file="solution/m"); write(0,*), m;
!    open(1,file="solution/f"); write(1,*), f;
    call proxJ(Proxm,Proxf,m,f,gamma,N,Q);
!    print *, "m : "
!    do i = 1,Q+1
!        print *, m(:,i), "ENDL";
!    end do
!    print *, "f : "
!    do i = 1,Q+1
!        print *, f(:,i), "ENDL";
!    end do
!    print *, "Proximite J m : "
!    do i = 1,Q+1
!        print *, Proxm(:,i), "ENDL";
!    end do
!    print *, "Proximite J f : "
!    do i = 1,Q+1
!        print *, Proxf(:,i), "ENDL";
!    end do
    
    !! Check prox G1 
    call proxG1(Umbar,Ufbar,Vm,Vf,mbar,fbar,m,f,gamma,N,Q);
!    print *, "Proximité G1 : Umbar :"
!    do i = 1,Q+1
!        print *, Umbar(:,i), "ENDL";
!    end do
!    print *, "Proximité G1 : Ufbar :"
!    do i = 1,Q+2
!        print *, Ufbar(:,i), "ENDL";
!    end do
!    print *, "Proximité G1 : Vm :"
!    do i = 1,Q+1
!        print *, Vm(:,i), "ENDL";
!    end do
!    print *, "Proximité G1 : Vf :"
!    do i = 1,Q+1
!        print *, Vf(:,i), "ENDL";
!    end do
    
    !! Check Projection on C 
    call projC(projCmbar,projCfbar,mbar,fbar,N,Q); 
!    print *, "Projection sur C mbar : ";
!    do i = 1,Q+1
!        print *, projCmbar(:,i), "ENDL";
!    end do 
!    print *, "Projection sur C fbar : ";
!    do i = 1,Q+2
!        print *, projCfbar(:,i), "ENDL";
!    end do
!    do while (1 .EQ. 1)
!    end do

    !! Check prox G2
!    mbar = 0; fbar = 0; m = 0; f= 0;
    call proxG2(mbart,fbart,mt,ft,mbar,fbar,m,f,N,Q);
!    print *, "Proximité G2 : mbart :"
!    do i = 1,Q+1
!        print *, mbart(:,i), "ENDL";
!    end do
!    print *, "Proximité G2 : fbart :"
!    do i = 1,Q+2
!        print *, fbart(:,i), "ENDL";
!    end do
!    print *, "Proximité G2 : mt :"
!    do i = 1,Q+1
!        print *, mt(:,i), "ENDL";
!    end do
!    print *, "Proximité G2 : ft :"
!    do i = 1,Q+1
!        print *, ft(:,i), "ENDL";
!    end do

    !!! Fin Zone de tests !!!

end subroutine
