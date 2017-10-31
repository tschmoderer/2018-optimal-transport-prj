!! etant donné mbar et fbar, 
!! calcul la projection sur l'ensemble C 
!! Appel la résolution du pbm de poisson

subroutine projC(projCmbar,projCfbar,mbar,fbar,N,Q)
    implicit none
    integer, intent(in) :: N,Q;
    double precision, dimension(N+2,Q+1), intent(in) :: mbar;
    double precision, dimension(N+1,Q+2), intent(in) :: fbar;
    double precision, dimension(N+2,Q+1), intent(out) :: projCmbar;
    double precision, dimension(N+1,Q+2), intent(out) :: projCfbar;

    !! Constantes dans la projection !!
    double precision, dimension(N+1,Q+1) :: Cst;
    double precision, dimension(N+2,Q+1) :: Cstmbar;
    double precision, dimension(N+1,Q+2) :: Cstfbar;

    !! Second memebre de la projection !!
    double precision, dimension(N+1,Q+1) :: d; ! le centre 
    double precision, dimension(Q+1) :: mleft,mright; ! les frontières 
    double precision, dimension(N+1) :: fup,fdown;

    !! Solution du pbm de poisson !!
    double precision, dimension(N+1,Q+1) :: solution;
    double precision, dimension(N+2,Q+1) :: solutionmbar;
    double precision, dimension(N+1,Q+2) :: solutionfbar;

    integer :: i;

    !! Initialisation !! 
    Cst = 0; Cstmbar = 0; Cstfbar = 0;
    d = 0; mleft = 0; mright = 0; fup = 0; fdown = 0;
    solution = 0; solutionmbar = 0; solutionfbar = 0;

    open (unit=3,file="files/Y");
    read (3,*), Cst;
    close (3);

    call divergence_adjoint(Cstmbar,Cstfbar,Cst,N,Q);

    call divergence(mbar,fbar,d,N,Q);
    call boundary(mbar,fbar,mleft,mright,fup,fdown,N,Q);

    open(unit=4,file='files/d'); write(4,*), d; close(4);
    open(unit=7,file="files/mL"); write(7,*), mleft; close(7);
    open(unit=8,file="files/mR"); write(8,*), mright; close(8);
    open(unit=9,file="files/fU"); write(9,*), fup; close(9);
    open(unit=10,file="files/fD"); write(10,*), fdown; close(10);
    
    call system('FreeFem++ -v 0 poisson_2d.pde');

    open(unit=3,file="files/solution"); read(3,*), solution; close (3)

    call divergence_adjoint(solutionmbar,solutionfbar,solution,N,Q);
  
    projCmbar = mbar - solutionmbar + Cstmbar;
    projCfbar = fbar - solutionfbar + Cstfbar;    
end subroutine