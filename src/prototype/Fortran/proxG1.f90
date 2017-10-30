subroutine proxG1(Umbar,Ufbar,Vm,Vf,mbar,fbar,m,f,g,N,Q)
    implicit none
    integer, intent(in) :: N,Q;
    double precision, dimension(N+1,Q+1) :: m,f; 
    double precision, dimension(N+2,Q+1) :: mbar;
    double precision, dimension(N+1,Q+2) :: fbar; 
    real, intent(in) :: g;

    !! La projection !!
    double precision, dimension(N+1,Q+1) :: Vm,Vf; ! projection de m sur J 
    double precision, dimension(N+2,Q+1) :: Umbar; ! projection de mbar sur C
    double precision, dimension(N+1,Q+2) :: Ufbar; ! projection de fbar sur C

    call projC(Umbar,Ufbar,mbar,fbar,N,Q);
    call proxJ(Vm,Vf,m,f,g,N,Q);
end subroutine