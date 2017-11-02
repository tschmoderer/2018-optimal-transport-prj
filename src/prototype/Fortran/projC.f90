!! etant donné mbar et fbar, 
!! calcul la projection sur l'ensemble C 
!! Appel la résolution du pbm de poisson

subroutine projC(projCmbar,projCfbar,mbar,fbar,N,Q)
    implicit none
    integer, intent(in) :: N,Q;
    double precision, dimension(Q+1,N+2), intent(in) :: mbar;
    double precision, dimension(Q+2,N+1), intent(in) :: fbar;
    double precision, dimension(Q+1,N+2), intent(out) :: projCmbar;
    double precision, dimension(Q+2,N+1), intent(out) :: projCfbar;

    !! Constantes dans la projection !!
    double precision, dimension(Q+1,N+1) :: Cst;
    double precision, dimension(Q+1,N+2) :: Cstmbar;
    double precision, dimension(Q+2,N+1) :: Cstfbar;

    !! Second memebre de la projection !!
    double precision, dimension(Q+1,N+1) :: d; ! le centre 
    double precision, dimension(Q+1) :: mleft,mright; ! les frontières 
    double precision, dimension(N+1) :: fup,fdown;

    !! Solution du pbm de poisson !!
    double precision, dimension(Q+1,N+1) :: solution;
    double precision, dimension(Q+1,N+2) :: solutionmbar;
    double precision, dimension(Q+2,N+1) :: solutionfbar;
    
    !test 
    double precision, dimension(0:N) :: GcX;
    integer :: i,j;

    !! Initialisation !! 

    Cst = 0; Cstmbar = 0; Cstfbar = 0;
    d = 0; mleft = 0; mright = 0; fup = 0; fdown = 0;
    solution = 0; solutionmbar = 0; solutionfbar = 0;

    open (unit=3,file="files/Y");
    do i = 1,Q+1
        do j = 1,N+1
            read(3,*), Cst(i,j);
        end do
    end do 
    close (3); 

    call divergence_adjoint(Cstmbar,Cstfbar,Cst,N,Q);

    call divergence(mbar,fbar,d,N,Q);
    call boundary(mbar,fbar,mleft,mright,fup,fdown,N,Q);
    ! test 
    do i=0,N
     GcX(i) = i/(1.0*N);
    end do
!    call finitial(GcX,N,fdown);
!    call ffinal(GcX,N,fup);

    open(unit=7,file="files/mL"); write(7,*), mleft; close(7);
    open(unit=8,file="files/mR"); write(8,*), mright; close(8);
    open(unit=9,file="files/fU"); write(9,*), fup; close(9);
    open(unit=10,file="files/fD"); write(10,*), fdown; close(10);
    
    open(unit=4,file='files/d'); 
    do i = 1,Q+1
        do j = 1,N+1
            write(4,*), d(i,j); 
        end do 
    end do 
    close(4);


    call system('FreeFem++ -v 0 poisson_2d.pde');

    open(unit=3,file="files/solution"); 
    do i = 1,Q+1
        do j = 1,N+1
            read(3,*), solution(i,j); 
        end do
    end do
    close (3)

    call divergence_adjoint(solutionmbar,solutionfbar,solution,N,Q);
  
    projCmbar = mbar - solutionmbar + Cstmbar;
    projCfbar = fbar - solutionfbar + Cstfbar;    

    !! de plus 

  !  projCmbar(1,:) = mbar(1,:) - solutionmbar(1,:) + Cstmbar(1,:);
  !  projCmbar(N+2,:) = mbar(N+2,:) - solutionmbar(N+2,:) + Cstmbar(N+2,:);
  !  projCfbar(:,1) = fbar(:,1) - solutionfbar(:,1) + Cstfbar(:,1);
  !  projCfbar(:,Q+2) = fbar(:,Q+2) - solutionfbar(:,Q+2) + Cstfbar(:,Q+2);
end subroutine