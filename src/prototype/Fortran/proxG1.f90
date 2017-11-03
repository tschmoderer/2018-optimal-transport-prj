!! opérateur de proximité G1
!! in : mbar, fbar, m, f
!! out : Umbar Ufbar, Vm,Vf

subroutine proxG1(Umbar,Ufbar,Vm,Vf,mbar,fbar,m,f,g,N,Q)
    implicit none
    integer, intent(in) :: N,Q;
    double precision, dimension(Q+1,N+1), intent(in) :: m,f; 
    double precision, dimension(Q+1,N+2), intent(in) :: mbar;
    double precision, dimension(Q+2,N+1), intent(in) :: fbar; 
    double precision, intent(in) :: g;

    !! La projection !!
    double precision, dimension(Q+1,N+1), intent(out) :: Vm,Vf; ! projection de m sur J 
    double precision, dimension(Q+1,N+2), intent(out) :: Umbar; ! projection de mbar sur C
    double precision, dimension(Q+2,N+1), intent(out) :: Ufbar; ! projection de fbar sur C

    call projC(Umbar,Ufbar,mbar,fbar,N,Q);
    call proxJ(Vm,Vf,m,f,g,N,Q);
end subroutine
