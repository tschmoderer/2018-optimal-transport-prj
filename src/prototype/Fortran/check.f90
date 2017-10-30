subroutine check(alpha,beta,gamma,N,Q)
    implicit none 
    double precision, intent(in) :: alpha,beta,gamma;
    integer, intent(in) :: N,Q;
    integer :: i,j;

    double precision, dimension(N+1,Q+1) :: m,f;
    double precision, dimension(N+2,Q+1) :: mbar;
    double precision, dimension(N+1,Q+2) :: fbar;

    !! Le cout de la solution !!
    double precision R;
    !! Divergence !! 
    double precision, dimension(N+1,Q+1) :: d;
    !! Boundary !!
    double precision, dimension(Q+1) :: mleft,mright;
    double precision, dimension(N+1) :: fup,fdown;
    !! Interpolsation !!
    double precision, dimension(N+1,Q+1) :: Interpm,Interpf;
    !! Prox J !!
    double precision, dimension(N+1,Q+1) :: Proxm,Proxf; 
    !! Proj C !!
    double precision, dimension(N+2,Q+1) :: projCmbar;
    double precision, dimension(N+1,Q+2) :: projCfbar;  
    !! Prox G1 !!
    double precision, dimension(N+1,Q+1) :: Vm,Vf;
    double precision, dimension(N+2,Q+1) :: Umbar;
    double precision, dimension(N+1,Q+2) :: Ufbar;
    !! Prox G2 !!
    double precision, dimension(N+1,Q+1) :: mt,ft;
    double precision, dimension(N+2,Q+1) :: mbart;
    double precision, dimension(N+1,Q+2) :: fbart;

    !!! Zone de tests !!!
    call RANDOM_NUMBER(m);
    call RANDOM_NUMBER(f);
    call RANDOM_NUMBER(mbar);
    call RANDOM_NUMBER(fbar);

    open(0,file="../m"); write(0,*) m; close(0);
    open(0,file="../f"); write(0,*) f; close(0);
    open(0,file="../mbar"); write(0,*) mbar; close(0);
    open(0,file="../fbar"); write(0,*) fbar; close(0);

    print *, "mbar : "
    do i = 1,Q+1
        print *, mbar(:,i), "ENDL";
    end do

    !! Check cost !!

    call cost(R,m,f,N,Q);
    print *, "Cout de la solution : ", R;

    !! Check divergence !!
    call divergence(mbar,fbar,d,N,Q);
    print *, "Divergence : ";
    do i = 1,Q+1
      print *, d(:,i), "ENDL";
    end do

    !! Check Boundary !!
    call boundary(mbar,fbar,mleft,mright,fup,fdown,N,Q)
    print *, "mbar : "
    do i = 1,Q+1
        print *, mbar(:,i), "ENDL";
    end do
    print *, "mleft : ", mleft; 
     print *, "fbar : "
    do i = 1,Q+2
        print *, fbar(:,i), "ENDL";
    end do
    print *, "fup : ", fup;

    !! Check Interpolation 
    call interpolation(mbar,fbar,Interpm,Interpf,N,Q)
    print *, "Interpolation : ";
    do i = 1,Q+1
        print *, Interpf(:,i), "ENDL";
    end do

    !! check adjoint
    call check_adjoint(N,Q)

    !! Check Proximity operator of J
    open(0,file="solution/m"); write(0,*), m;
    open(1,file="solution/f"); write(1,*), f;
    call proxJ(Proxm,Proxf,m,f,gamma,N,Q);

    !! Check Projection on C 
    call projC(projCmbar,projCfbar,mbar,fbar,N,Q); 

    !! Check prox G1 
    call proxG1(Umbar,Ufbar,Vm,Vf,mbar,fbar,m,f,gamma,N,Q);

    !! Check prox G2
    call proxG2(mbart,fbart,mt,ft,mbar,fbar,m,f,N,Q);

    !!! Fin Zone de tests !!!

end subroutine