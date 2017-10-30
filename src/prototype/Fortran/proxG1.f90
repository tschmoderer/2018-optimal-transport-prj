subroutine proxG1(Umbar,Ufbar,Vm,Vf,mbar,fbar,m,f,g,N,Q)
    implicit none
    integer, intent(in) :: N,Q;
    double precision, dimension(N+1,Q+1), intent(in) :: m,f; 
    double precision, dimension(N+2,Q+1), intent(in) :: mbar;
    double precision, dimension(N+1,Q+2), intent(in) :: fbar; 
    real, intent(in) :: g;
    
    !! La projection !!
    double precision, dimension(N+1,Q+1), intent(out) :: Vm,Vf; ! projection de m sur J 
    double precision, dimension(N+2,Q+1), intent(out) :: Umbar; ! projection de mbar sur C
    double precision, dimension(N+1,Q+2), intent(out) :: Ufbar; ! projection de fbar sur C

    call projC(Umbar,Ufbar,mbar,fbar,N,Q);
    call proxJ(Vm,Vf,m,f,g,N,Q);
end subroutine