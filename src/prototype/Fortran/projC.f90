!! etant donné mbar et fbar, 
!! calcul la projection sur l'ensemble C 
!! Appel la résolution du pbm de poisson

subroutine projC(projCmbar,projCfbar,mbar,fbar,N,Q)
    implicit none
    integer, intent(in) :: N,Q;
    real, dimension(N+2,Q+1), intent(in) :: mbar;
    real, dimension(N+1,Q+2), intent(in) :: fbar;
    real, dimension(N+2,Q+1), intent(out) :: projCmbar;
    real, dimension(N+1,Q+2), intent(out) :: projCfbar;

    !! Constantes dans la projection !!
    real, dimension(N+1,Q+1) :: Cst;
    real, dimension(N+2,Q+1) :: Cstmbar;
    real, dimension(N+1,Q+2) :: Cstfbar;

    !! Second memebre de la projection !!
    real, dimension(N+1,Q+1) :: d; ! le centre 
    real, dimension(Q+1) :: mleft,mright;
    real, dimension(N+1) :: fup,fdown;

    !! Solution du pbm de poisson !!
    real, dimension(N+1,Q+1) :: solution;
    real, dimension(N+2,Q+1) :: solutionmbar;
    real, dimension(N+1,Q+2) :: solutionfbar;

    !! Initialisation !! 
    Cst = 0; Cstmbar = 0; Cstfbar = 0;
    d = 0; mleft = 0; mright = 0; fup = 0; fdown = 0;
    solution = 0; solutionmbar = 0; solutionfbar = 0;

    open (unit=3,file="files/Y");
    read (3,*), Cst;
    close (unit=3);

    call divergence_adjoint(Cstmbar,Cstfbar,Cst,N,Q);

    call divergence(mbar,fbar,d,N,Q);
    call boundary(mbar,fbar,mleft,mright,fup,fdown,N,Q);

    open(unit=4,file='files/d'); 
    open(unit=7,file="files/mL");
    open(unit=8,file="files/mR");
    open(unit=9,file="files/fU");
    open(unit=10,file="files/fD");
    
    write(4,*), d;
    write(7,*), mleft;
    write(8,*), mright;
    write(9,*), fup;
    write(10,*), fdown;
    
    call system('FreeFem++ poisson_2d.pde');

    open(11,file="files/solution");
    
    read(11,*), solution;
 
    call divergence_adjoint(solutionmbar,solutionmbar,solution,N,Q);
  

  !  projCmbar = mbar - solutionmbar + Cstmbar;
  !  projCfbar = fbar - solutionfbar + Cstfbar;
    projCmbar = 0;
    projCfbar = 0;

  
end subroutine